# Pressure advance calibration from load cell force
#
# Sweeps pressure advance and picks the last value before the force
# trace undershoots after a fast->slow flow step.
#
# Pressure advance is only applied to moves that also move X or Y, see
# PrinterExtruder.process_move. So this extrudes while travelling at a
# constant XY speed and steps the extrusion width rather than the feed
# rate. Constant XY speed keeps the toolhead's own inertial load out of
# the measurement; only the extrusion width changes across the step.
#
# This file may be distributed under the terms of the GNU GPLv3 license.
from __future__ import annotations

import csv
import os
import statistics
import time

from klippy.extras.load_cell.load_cell import SampleStructure
from klippy.extras.load_cell_flow_test import find_load_cell

# Force is smoothed over this long before peaks are taken, because the
# minimum of a noisy window is biased low by the noise itself.
SMOOTH_SECONDS = 0.01
# The flow step must be this many times the idle noise to be usable.
MIN_STEP_SIGMA = 20.0
# Fraction of each segment treated as settled, measured from its end.
SETTLED_FRACTION = 0.3
# Fraction of the trailing segment searched for the undershoot dip.
DIP_FRACTION = 0.6


def boxcar(values, width):
    if width < 3:
        return list(values)
    half = width // 2
    prefix = [0.0]
    for value in values:
        prefix.append(prefix[-1] + value)
    out = []
    for idx in range(len(values)):
        lo = max(0, idx - half)
        hi = min(len(values), idx + half + 1)
        out.append((prefix[hi] - prefix[lo]) / (hi - lo))
    return out


class Trace:
    """Smoothed, sign-normalised force over one sweep pass."""

    def __init__(self, samples, sign, width):
        column = SampleStructure.COUNTS.value
        self.times = [s[SampleStructure.TIME.value] for s in samples]
        raw = [sign * s[column] for s in samples]
        self.values = boxcar(raw, width)

    def window(self, start, end):
        return [
            v
            for t, v in zip(self.times, self.values)
            if start <= t <= end
        ]

    def settled(self, start, end):
        span = end - start
        values = self.window(end - span * SETTLED_FRACTION, end)
        return statistics.fmean(values) if values else None

    def dip(self, start, end):
        span = end - start
        values = self.window(start, start + span * DIP_FRACTION)
        return min(values) if values else None

    def peak(self, start, end):
        span = end - start
        values = self.window(start, start + span * DIP_FRACTION)
        return max(values) if values else None


class PassResult:
    def __init__(self, advance):
        self.advance = advance
        self.step = 0.0
        self.undershoot_pct = 0.0
        self.overshoot_pct = 0.0


class PressureAdvanceCalibrate:
    def __init__(self, config):
        self.printer = config.get_printer()
        gcode = self.printer.lookup_object("gcode")
        gcode.register_command(
            "CALIBRATE_PRESSURE_ADVANCE",
            self.cmd_CALIBRATE_PRESSURE_ADVANCE,
            desc=self.cmd_CALIBRATE_PRESSURE_ADVANCE_help,
        )

    cmd_CALIBRATE_PRESSURE_ADVANCE_help = (
        "Sweep pressure advance and pick a value from load cell force"
    )

    def cmd_CALIBRATE_PRESSURE_ADVANCE(self, gcmd):
        lc_name, load_cell = find_load_cell(self.printer, gcmd)
        toolhead = self.printer.lookup_object("toolhead")
        extruder = toolhead.get_extruder()
        if not extruder.get_name():
            raise gcmd.error("No extruder configured")
        status = toolhead.get_status(self.printer.get_reactor().monotonic())
        if not all(a in status.get("homed_axes", "") for a in "xyz"):
            raise gcmd.error("Home all axes first")

        start = gcmd.get_float("START", 0.0, minval=0.0)
        stop = gcmd.get_float("STOP", 0.1, above=0.0)
        step = gcmd.get_float("STEP", 0.005, above=0.0)
        if stop <= start:
            raise gcmd.error("STOP must be greater than START")
        flow_low = gcmd.get_float("FLOW_LOW", 4.0, above=0.0)
        flow_high = gcmd.get_float("FLOW_HIGH", 6.0, above=0.0)
        if flow_high <= flow_low:
            raise gcmd.error("FLOW_HIGH must be greater than FLOW_LOW")
        xy_speed = gcmd.get_float("TRAVEL_SPEED", 5.0, above=0.0)
        duration = gcmd.get_float("DURATION", 0.8, minval=0.2, maxval=5.0)
        settle = gcmd.get_float("SETTLE", 0.5, minval=0.2, maxval=5.0)
        lead = gcmd.get_float("LEAD", 0.3, minval=0.1, maxval=5.0)
        threshold = gcmd.get_float("THRESHOLD", 2.0, above=0.0)

        area = extruder.filament_area
        cross_high = flow_high / xy_speed
        # Slow travel keeps the purge in one place, which means a fat
        # extrusion. max_extrude_cross_section exists to catch slicer
        # mistakes and is not meaningful here, so raise it for the sweep
        # and put it back afterwards.
        old_nozzle = extruder.nozzle_diameter
        old_cross = extruder.max_extrude_ratio * area
        need_cross = cross_high * 1.5
        # The lookahead limits junction speed between the flow segments
        # to instantaneous_corner_velocity / the change in extrusion
        # ratio. Travel speed cancels out of that, leaving a condition on
        # the flow difference alone. Breach it and XY decelerates right
        # at the step, putting an inertial transient on the edge being
        # measured and smearing the step itself.
        max_flow_step = extruder.instant_corner_v * area
        if flow_high - flow_low > max_flow_step:
            raise gcmd.error(
                "A %.2f mm^3/s flow step makes the lookahead slow XY down at "
                "the transition, which smears the edge and adds an inertial "
                "transient to it. Keep FLOW_HIGH - FLOW_LOW under %.2f "
                "mm^3/s, or raise instantaneous_corner_velocity in "
                "[extruder]." % (flow_high - flow_low, max_flow_step)
            )

        # G1 targets are in gcode coordinates, the axis limits are in
        # toolhead coordinates, and a G92 offset separates the two.
        gcode_move = self.printer.lookup_object("gcode_move")
        gc_pos = gcode_move.get_status()["gcode_position"]
        th_pos = toolhead.get_position()
        x0 = gcmd.get_float("X", gc_pos[0])
        y0 = gcmd.get_float("Y", gc_pos[1])
        z0 = gcmd.get_float("Z", gc_pos[2])
        seg_len = xy_speed * duration
        lead_len = xy_speed * lead
        x_offset = th_pos[0] - gc_pos[0]
        x_end = x0 + 3.0 * seg_len + 2.0 * lead_len
        x_max = status.get("axis_maximum")
        if x_max is not None and x_end + x_offset > x_max[0]:
            raise gcmd.error(
                "Sweep needs X from %.1f to %.1f, past the %.1f limit. Lower "
                "DURATION or TRAVEL_SPEED, or pass a smaller X."
                % (x0 + x_offset, x_end + x_offset, x_max[0])
            )

        advances = []
        value = start
        while value <= stop + 1e-9:
            advances.append(round(value, 6))
            value += step
        gcmd.respond_info(
            "Sweeping %i pressure advance values on '%s': %.4f to %.4f. "
            "%.1f->%.1f mm^3/s at %.0fmm/s, %.1fmm per segment plus %.1fmm "
            "lead in/out, X %.1f to %.1f, Z=%.1f. Extruding into free air."
            % (
                len(advances),
                lc_name,
                advances[0],
                advances[-1],
                flow_low,
                flow_high,
                xy_speed,
                seg_len,
                lead_len,
                x0,
                x_end,
                z0,
            )
        )

        gcode = self.printer.lookup_object("gcode")
        gcode.run_script_from_command("SAVE_GCODE_STATE NAME=_pa_calibrate")
        results = []
        rows = []
        try:
            if need_cross > old_cross:
                gcode.run_script_from_command(
                    "SET_NOZZLE_DIAMETER DIAMETER=%.4f "
                    "MAX_EXTRUDE_CROSS_SECTION=%.4f"
                    % (old_nozzle, need_cross)
                )
            gcode.run_script_from_command("G90")
            gcode.run_script_from_command("M83")
            for advance in advances:
                result, pass_rows = self._run_pass(
                    gcmd,
                    load_cell,
                    toolhead,
                    gcode,
                    advance,
                    x0,
                    y0,
                    z0,
                    seg_len,
                    lead_len,
                    xy_speed,
                    area,
                    flow_low,
                    flow_high,
                    duration,
                    settle,
                    lead,
                )
                results.append(result)
                rows.extend(pass_rows)
        finally:
            if need_cross > old_cross:
                gcode.run_script_from_command(
                    "SET_NOZZLE_DIAMETER DIAMETER=%.4f "
                    "MAX_EXTRUDE_CROSS_SECTION=%.4f"
                    % (old_nozzle, old_cross)
                )
            gcode.run_script_from_command(
                "RESTORE_GCODE_STATE NAME=_pa_calibrate MOVE=0"
            )
        self._write_csv(gcmd, rows)
        self._report(gcmd, results, threshold)

    def _write_csv(self, gcmd, rows):
        log_file = self.printer.get_start_args().get("log_file")
        directory = os.path.dirname(log_file) if log_file else "/tmp"
        name = "pa_calibrate_%s.csv" % (time.strftime("%Y%m%d_%H%M%S"),)
        path = os.path.join(directory, name)
        try:
            with open(path, "w", newline="") as handle:
                writer = csv.writer(handle)
                writer.writerow(
                    ["advance", "time", "counts", "force_g", "segment"]
                )
                writer.writerows(rows)
        except OSError as err:
            gcmd.respond_info("Could not write %s: %s" % (path, err))
            return
        gcmd.respond_info("Raw force trace written to %s" % (path,))

    def _run_pass(
        self,
        gcmd,
        load_cell,
        toolhead,
        gcode,
        advance,
        x0,
        y0,
        z0,
        seg_len,
        lead_len,
        xy_speed,
        area,
        flow_low,
        flow_high,
        duration,
        settle,
        lead,
    ):
        gcode.run_script_from_command(
            "SET_PRESSURE_ADVANCE ADVANCE=%.6f" % (advance,)
        )
        gcode.run_script_from_command(
            "G1 X%.3f Y%.3f Z%.3f F%.0f" % (x0, y0, z0, xy_speed * 60.0)
        )
        feed = xy_speed * 60.0
        # Lead in and lead out at the same speed and the low flow, so
        # that acceleration and deceleration happen outside every
        # measured window. Without them the settled window of the last
        # segment sits in the decel ramp and reads toolhead inertia
        # instead of melt pressure.
        plan = (
            ("lead_end", lead_len, flow_low * lead / area),
            ("slow_end", seg_len, flow_low * duration / area),
            ("fast_end", seg_len, flow_high * duration / area),
            ("slow2_end", seg_len, flow_low * duration / area),
            ("out_end", lead_len, flow_low * lead / area),
        )

        collector = load_cell.get_collector()
        marks = {}

        def mark(label):
            def callback(print_time):
                marks[label] = print_time

            toolhead.register_lookahead_callback(callback)

        t_start = toolhead.get_last_move_time()
        collector.start_collecting(min_time=t_start)
        # dwell() flushes, so the idle baseline is a clean stop
        toolhead.dwell(settle)
        marks["idle_end"] = toolhead.get_last_move_time()
        # The three extruding moves are queued without flushing between
        # them, so the lookahead blends them and the flow step stays
        # sharp. Lookahead callbacks give the exact boundary print time
        # without forcing a stop.
        x_target = x0
        for label, dist, e_dist in plan:
            x_target += dist
            gcode.run_script_from_command(
                "G1 X%.4f E%.4f F%.0f" % (x_target, e_dist, feed)
            )
            mark(label)
        toolhead.dwell(settle)
        t_end = toolhead.get_last_move_time()
        samples, errors = collector.collect_until(t_end)
        if not samples:
            raise gcmd.error("No samples returned from sensor!")
        load_cell.validate_samples(samples, errors)
        result = self._analyse(
            gcmd, load_cell, advance, samples, marks, t_start
        )
        return result, self._rows(advance, samples, marks, t_start)

    @staticmethod
    def _rows(advance, samples, marks, t_start):
        # Label each sample with the segment it falls in, so a plot can
        # show where the analysis windows actually sat.
        bounds = [
            ("idle", marks["idle_end"]),
            ("lead", marks["lead_end"]),
            ("slow", marks["slow_end"]),
            ("fast", marks["fast_end"]),
            ("slow2", marks["slow2_end"]),
            ("out", marks["out_end"]),
        ]
        rows = []
        for sample in samples:
            when = sample[SampleStructure.TIME.value]
            label = "tail"
            for name, end in bounds:
                if when <= end:
                    label = name
                    break
            rows.append(
                [
                    "%.6f" % (advance,),
                    "%.6f" % (when - t_start,),
                    sample[SampleStructure.COUNTS.value],
                    "%.3f" % (sample[SampleStructure.FORCE_G.value],),
                    label,
                ]
            )
        return rows

    def _analyse(self, gcmd, load_cell, advance, samples, marks, t_start):
        for label in (
            "idle_end",
            "lead_end",
            "slow_end",
            "fast_end",
            "slow2_end",
            "out_end",
        ):
            if label not in marks:
                raise gcmd.error(
                    "Lookahead callback '%s' never fired" % (label,)
                )
        column = SampleStructure.COUNTS.value
        idle = [
            s[column]
            for s in samples
            if t_start <= s[SampleStructure.TIME.value] <= marks["idle_end"]
        ]
        slow_raw = [
            s[column]
            for s in samples
            if marks["lead_end"]
            <= s[SampleStructure.TIME.value]
            <= marks["slow_end"]
        ]
        fast_raw = [
            s[column]
            for s in samples
            if marks["slow_end"]
            <= s[SampleStructure.TIME.value]
            <= marks["fast_end"]
        ]
        if len(idle) < 2 or not slow_raw or not fast_raw:
            raise gcmd.error("Too few samples in the sweep window")
        noise = statistics.pstdev(idle)
        sign = (
            1.0
            if statistics.fmean(fast_raw) > statistics.fmean(slow_raw)
            else -1.0
        )
        sps = load_cell.sensor.get_samples_per_second()
        trace = Trace(samples, sign, max(3, int(SMOOTH_SECONDS * sps)))

        fast_level = trace.settled(marks["slow_end"], marks["fast_end"])
        slow_level = trace.settled(marks["fast_end"], marks["slow2_end"])
        dip = trace.dip(marks["fast_end"], marks["slow2_end"])
        peak = trace.peak(marks["slow_end"], marks["fast_end"])
        if None in (fast_level, slow_level, dip, peak):
            raise gcmd.error("Empty analysis window - raise DURATION")
        step = fast_level - slow_level
        if noise > 0.0 and step < MIN_STEP_SIGMA * noise:
            raise gcmd.error(
                "Flow step is only %.1f sigma above noise (fast settled "
                "%.0f, slow settled %.0f, noise %.0f counts). A negative "
                "step means something other than flow is loading the cell "
                "in the measured windows."
                % (step / noise, fast_level, slow_level, noise)
            )
        result = PassResult(advance)
        result.step = step
        result.undershoot_pct = 100.0 * (slow_level - dip) / step
        result.overshoot_pct = 100.0 * (peak - fast_level) / step
        return result

    def _report(self, gcmd, results, threshold):
        gcmd.respond_info("advance   undershoot%   overshoot%   step(counts)")
        for res in results:
            gcmd.respond_info(
                "%7.4f   %10.2f   %9.2f   %12.0f"
                % (
                    res.advance,
                    res.undershoot_pct,
                    res.overshoot_pct,
                    res.step,
                )
            )
        over = [r for r in results if r.undershoot_pct >= threshold]
        if not over:
            gcmd.respond_info(
                "No undershoot above %.1f%% anywhere in the sweep. The right "
                "value is above STOP=%.4f - extend the range."
                % (threshold, results[-1].advance)
            )
            return
        first = over[0]
        if first is results[0]:
            gcmd.respond_info(
                "Already undershooting %.2f%% at the lowest value tested. "
                "The right value is below START=%.4f - extend the range "
                "downward." % (first.undershoot_pct, results[0].advance)
            )
            return
        chosen = results[results.index(first) - 1]
        gcmd.respond_info(
            "Undershoot crosses %.1f%% at advance %.4f, so the last good "
            "value is %.4f." % (threshold, first.advance, chosen.advance)
        )
        gcmd.respond_info(
            "Set it with: SET_PRESSURE_ADVANCE ADVANCE=%.4f" % (chosen.advance,)
        )


def load_config(config):
    return PressureAdvanceCalibrate(config)
