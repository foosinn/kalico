# Pressure advance calibration from load cell force
#
# Melt pressure takes seconds to equilibrate, far longer than one flow
# step, so a sweep that starts and stops extrusion per K value measures
# its own startup transient instead of pressure advance. This drives the
# melt into a periodic steady state instead: prime once, then cycle
# slow/fast continuously, never stopping, and read only the last cycle
# of each K.
#
# Pressure advance is only applied to moves that also move X or Y, see
# PrinterExtruder.process_move, so every segment travels. The travel
# speed is constant throughout, so the toolhead's own inertia is a
# constant offset; only the extrusion width steps. Travel direction
# alternates per K value to keep the purge in a bounded band.
#
# This file may be distributed under the terms of the GNU GPLv3 license.
from __future__ import annotations

import csv
import os
import random
import statistics
import time

from klippy.extras.load_cell.load_cell import SampleStructure
from klippy.extras.load_cell_flow_test import find_load_cell

# Force is smoothed over this long before peaks are taken, because the
# minimum of a noisy window is biased low by the noise itself.
SMOOTH_SECONDS = 0.01
# The flow step must be this many times the idle noise to be usable.
MIN_STEP_SIGMA = 20.0
# Fraction of a segment treated as settled, measured from its end.
SETTLED_FRACTION = 0.3
# Fraction of the trailing segment searched for the undershoot dip.
DIP_FRACTION = 0.6
# A settled window drifting by more than this fraction of the step over
# its own length has not settled, and its levels mean nothing.
MAX_TAIL_DRIFT = 0.1


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


def slope(points):
    """Least squares slope of (time, value) pairs."""
    if len(points) < 2:
        return 0.0
    mean_t = statistics.fmean([p[0] for p in points])
    mean_v = statistics.fmean([p[1] for p in points])
    num = sum((p[0] - mean_t) * (p[1] - mean_v) for p in points)
    den = sum((p[0] - mean_t) ** 2 for p in points)
    return num / den if den else 0.0


class Move:
    """One queued segment of the sweep."""

    def __init__(self, label, group, seconds, flow, advance):
        self.label = label
        self.group = group
        self.seconds = seconds
        self.flow = flow
        self.advance = advance
        self.start = 0.0
        self.end = 0.0


class PassResult:
    def __init__(self, advance, order):
        self.advance = advance
        self.order = order
        self.step = 0.0
        self.undershoot_pct = 0.0
        self.overshoot_pct = 0.0
        self.drift_pct = 0.0
        self.spread = 0.0
        self.edges = 0


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
        reactor = self.printer.get_reactor()
        status = toolhead.get_status(reactor.monotonic())
        if not all(a in status.get("homed_axes", "") for a in "xyz"):
            raise gcmd.error("Home all axes first")

        start = gcmd.get_float("START", 0.0, minval=0.0)
        stop = gcmd.get_float("STOP", 0.08, above=0.0)
        step = gcmd.get_float("STEP", 0.01, above=0.0)
        if stop <= start:
            raise gcmd.error("STOP must be greater than START")
        flow_low = gcmd.get_float("FLOW_LOW", 5.0, above=0.0)
        flow_high = gcmd.get_float("FLOW_HIGH", 30.0, above=0.0)
        if flow_high <= flow_low:
            raise gcmd.error("FLOW_HIGH must be greater than FLOW_LOW")
        xy_speed = gcmd.get_float("TRAVEL_SPEED", 2.0, above=0.0)
        duration = gcmd.get_float("DURATION", 2.0, minval=0.2, maxval=10.0)
        cycles = gcmd.get_int("CYCLES", 5, minval=2, maxval=20)
        prime = gcmd.get_float("PRIME", 30.0, minval=0.0, maxval=600.0)
        # Alternating keeps the purge in one band, but a few grams of
        # free-air extrudate coils into a tangle tall enough to load the
        # nozzle and trip the probe's force_safety_limit. Spreading it
        # along one line costs bed width and avoids that.
        alternate = gcmd.get_int("ALTERNATE", 0, minval=0, maxval=1)
        settle = gcmd.get_float("SETTLE", 0.5, minval=0.2, maxval=5.0)
        threshold = gcmd.get_float("THRESHOLD", 3.0, above=0.0)
        seed = gcmd.get_int("SEED", 0)

        area = extruder.filament_area
        # The lookahead caps junction speed between flow segments at
        # instantaneous_corner_velocity / the change in extrusion ratio.
        # Travel speed cancels out, leaving a limit on the flow step.
        # Pressure advance injects pa * (change in extruder velocity), so
        # the flow step is what sets the signal, not the absolute flow.
        # instantaneous_corner_velocity caps that step: breach it and the
        # lookahead decelerates XY at the transition, smearing the edge
        # and putting an inertial transient on it. It is only read by
        # PrinterExtruder.calc_junction, so raise it for the sweep and
        # put it back afterwards rather than making the user edit config.
        old_corner_v = extruder.instant_corner_v
        need_corner_v = (flow_high - flow_low) / area * 1.2
        if need_corner_v > old_corner_v:
            gcmd.respond_info(
                "Raising instantaneous_corner_velocity %.2f -> %.2f mm/s for "
                "the sweep so a %.1f mm^3/s flow step does not slow XY at "
                "the transition."
                % (old_corner_v, need_corner_v, flow_high - flow_low)
            )
        # Slow travel keeps the purge local, which means a fat extrusion.
        # max_extrude_cross_section guards against slicer mistakes and is
        # not meaningful here, so raise it and put it back afterwards.
        old_nozzle = extruder.nozzle_diameter
        old_cross = extruder.max_extrude_ratio * area
        need_cross = (flow_high / xy_speed) * 1.5

        advances = []
        value = start
        while value <= stop + 1e-9:
            advances.append(round(value, 6))
            value += step
        # K correlated with time is indistinguishable from slow drift, so
        # run them in a shuffled but reproducible order.
        order = list(advances)
        random.Random(seed).shuffle(order)

        gcode_move = self.printer.lookup_object("gcode_move")
        gc_pos = gcode_move.get_status()["gcode_position"]
        th_pos = toolhead.get_position()
        x0 = gcmd.get_float("X", gc_pos[0])
        y0 = gcmd.get_float("Y", gc_pos[1])
        z0 = gcmd.get_float("Z", gc_pos[2])

        plan = self._build_plan(order, prime, cycles, duration, flow_low,
                                flow_high)
        self._check_bounds(gcmd, plan, x0, xy_speed,
                           th_pos[0] - gc_pos[0], status.get("axis_minimum"),
                           status.get("axis_maximum"), alternate)
        total = sum(m.seconds for m in plan)
        gcmd.respond_info(
            "Sweeping %i pressure advance values on '%s': %.4f to %.4f in "
            "shuffled order, %i cycles each. %.1f->%.1f mm^3/s at %.1fmm/s, "
            "%.0fs of continuous extrusion, Z=%.1f. Extruding into free air."
            % (len(advances), lc_name, advances[0], advances[-1], cycles,
               flow_low, flow_high, xy_speed, total + prime, z0)
        )

        gcode = self.printer.lookup_object("gcode")
        gcode.run_script_from_command("SAVE_GCODE_STATE NAME=_pa_calibrate")
        try:
            if need_cross > old_cross:
                gcode.run_script_from_command(
                    "SET_NOZZLE_DIAMETER DIAMETER=%.4f "
                    "MAX_EXTRUDE_CROSS_SECTION=%.4f" % (old_nozzle, need_cross)
                )
            if need_corner_v > old_corner_v:
                extruder.instant_corner_v = need_corner_v
            gcode.run_script_from_command("G90")
            gcode.run_script_from_command("M83")
            samples, idle_end, t_start = self._run(
                gcmd, load_cell, toolhead, gcode, plan, x0, y0, z0,
                xy_speed, area, settle, alternate
            )
        finally:
            extruder.instant_corner_v = old_corner_v
            if need_cross > old_cross:
                gcode.run_script_from_command(
                    "SET_NOZZLE_DIAMETER DIAMETER=%.4f "
                    "MAX_EXTRUDE_CROSS_SECTION=%.4f" % (old_nozzle, old_cross)
                )
            gcode.run_script_from_command(
                "RESTORE_GCODE_STATE NAME=_pa_calibrate MOVE=0"
            )
        self._write_csv(gcmd, plan, samples, t_start)
        results = self._analyse(
            gcmd, load_cell, plan, samples, idle_end, t_start, order, cycles
        )
        self._report(gcmd, results, threshold)

    @staticmethod
    def _build_plan(order, prime, cycles, duration, flow_low, flow_high):
        plan = []
        # Prime with the same slow/fast cycle that will be measured, not
        # a steady low flow. The melt equilibrates to the mean flow, so
        # priming at flow_low leaves the first K values climbing toward
        # the cycle's mean and that climb can outrun the step itself.
        # One group per cycle keeps the travel in a band.
        remaining = prime
        index = 0
        while remaining > 1e-6:
            for kind, flow in (("slow", flow_low), ("fast", flow_high)):
                if remaining <= 1e-6:
                    break
                seconds = min(duration, remaining)
                plan.append(
                    Move(
                        "prime%d_%s" % (index, kind),
                        "p%d" % (index,),
                        seconds,
                        flow,
                        order[0],
                    )
                )
                remaining -= seconds
            index += 1
        for k_idx, advance in enumerate(order):
            for cycle in range(cycles):
                for kind, flow in (("slow", flow_low), ("fast", flow_high)):
                    plan.append(
                        Move(
                            "k%d_c%d_%s" % (k_idx, cycle, kind),
                            "k%d" % (k_idx,),
                            duration,
                            flow,
                            advance,
                        )
                    )
        return plan

    @staticmethod
    def _check_bounds(gcmd, plan, x0, xy_speed, offset, axis_min, axis_max,
                      alternate):
        # Direction flips per K value, so walk the plan and find the
        # extremes rather than assuming a single straight line.
        x = x0
        lo = hi = x0
        direction = 1.0
        current = None
        for move in plan:
            if alternate and current is not None and move.group != current:
                direction = -direction
            current = move.group
            x += direction * xy_speed * move.seconds
            lo = min(lo, x)
            hi = max(hi, x)
        if axis_max is not None and hi + offset > axis_max[0]:
            raise gcmd.error(
                "Sweep reaches X=%.1f, past the %.1f limit. Lower CYCLES, "
                "DURATION or TRAVEL_SPEED, or pass a smaller X."
                % (hi + offset, axis_max[0])
            )
        if axis_min is not None and lo + offset < axis_min[0]:
            raise gcmd.error(
                "Sweep reaches X=%.1f, below the %.1f limit. Lower CYCLES, "
                "DURATION or TRAVEL_SPEED, or pass a larger X."
                % (lo + offset, axis_min[0])
            )

    def _run(self, gcmd, load_cell, toolhead, gcode, plan, x0, y0, z0,
             xy_speed, area, settle, alternate):
        feed = xy_speed * 60.0
        gcode.run_script_from_command(
            "G1 X%.3f Y%.3f Z%.3f F%.0f" % (x0, y0, z0, feed)
        )
        collector = load_cell.get_collector()
        t_start = toolhead.get_last_move_time()
        collector.start_collecting(min_time=t_start)
        toolhead.dwell(settle)
        idle_end = toolhead.get_last_move_time()

        def mark(move):
            def callback(print_time):
                move.end = print_time

            toolhead.register_lookahead_callback(callback)

        x = x0
        direction = 1.0
        current = None
        previous_advance = None
        for move in plan:
            if alternate and current is not None and move.group != current:
                direction = -direction
            current = move.group
            if move.advance != previous_advance:
                # Goes through a lookahead callback unless smooth_time
                # changes, so motion keeps running across the change.
                gcode.run_script_from_command(
                    "SET_PRESSURE_ADVANCE ADVANCE=%.6f" % (move.advance,)
                )
                previous_advance = move.advance
            x += direction * xy_speed * move.seconds
            e_dist = move.flow * move.seconds / area
            gcode.run_script_from_command(
                "G1 X%.4f E%.4f F%.0f" % (x, e_dist, feed)
            )
            mark(move)
        toolhead.dwell(settle)
        t_end = toolhead.get_last_move_time()
        samples, errors = collector.collect_until(t_end)
        if not samples:
            raise gcmd.error("No samples returned from sensor!")
        load_cell.validate_samples(samples, errors)
        # start times are only known once every callback has fired
        previous = idle_end
        for move in plan:
            if not move.end:
                raise gcmd.error(
                    "Lookahead callback for '%s' never fired" % (move.label,)
                )
            move.start = previous
            previous = move.end
        return samples, idle_end, t_start

    def _analyse(self, gcmd, load_cell, plan, samples, idle_end, t_start,
                 order, cycles):
        column = SampleStructure.COUNTS.value
        idle = [
            s[column]
            for s in samples
            if t_start <= s[SampleStructure.TIME.value] <= idle_end
        ]
        if len(idle) < 2:
            raise gcmd.error("Too few idle samples for a noise estimate")
        noise = statistics.pstdev(idle)
        by_label = {m.label: m for m in plan}
        sps = load_cell.sensor.get_samples_per_second()
        times = [s[SampleStructure.TIME.value] for s in samples]
        width = max(3, int(SMOOTH_SECONDS * sps))

        probe = by_label["k0_c%d_fast" % (cycles - 2,)]
        sign = 1.0
        window = [
            s[column]
            for s in samples
            if probe.start <= s[SampleStructure.TIME.value] <= probe.end
        ]
        slow_probe = by_label["k0_c%d_slow" % (cycles - 1,)]
        slow_window = [
            s[column]
            for s in samples
            if slow_probe.start <= s[SampleStructure.TIME.value]
            <= slow_probe.end
        ]
        if window and slow_window:
            if statistics.fmean(window) < statistics.fmean(slow_window):
                sign = -1.0
        values = boxcar([sign * s[column] for s in samples], width)
        pairs = list(zip(times, values))

        def window_of(move, lo=0.0, hi=1.0):
            span = move.end - move.start
            a = move.start + span * lo
            b = move.start + span * hi
            return [p for p in pairs if a <= p[0] <= b]

        # Skip the cycles right after a K change, then average every
        # remaining edge. Scatter falls roughly as sqrt(edges).
        skip = min(2, cycles - 2)
        results = []
        for k_idx, advance in enumerate(order):
            steps, unders, overs, drifts = [], [], [], []
            for cycle in range(skip, cycles - 1):
                fast = by_label["k%d_c%d_fast" % (k_idx, cycle)]
                slow = by_label["k%d_c%d_slow" % (k_idx, cycle + 1)]
                fast_tail = window_of(fast, 1.0 - SETTLED_FRACTION, 1.0)
                slow_tail = window_of(slow, 1.0 - SETTLED_FRACTION, 1.0)
                dip_window = window_of(slow, 0.0, DIP_FRACTION)
                peak_window = window_of(fast, 0.0, DIP_FRACTION)
                if not (fast_tail and slow_tail and dip_window
                        and peak_window):
                    raise gcmd.error(
                        "Empty analysis window at advance %.4f - raise "
                        "DURATION" % (advance,)
                    )
                fast_level = statistics.fmean([p[1] for p in fast_tail])
                slow_level = statistics.fmean([p[1] for p in slow_tail])
                step = fast_level - slow_level
                if noise > 0.0 and step < MIN_STEP_SIGMA * noise:
                    raise gcmd.error(
                        "Flow step at advance %.4f is only %.1f sigma above "
                        "noise (fast %.0f, slow %.0f, noise %.0f counts)."
                        % (advance, step / noise, fast_level, slow_level,
                           noise)
                    )
                steps.append(step)
                unders.append(
                    100.0 * (slow_level - min(p[1] for p in dip_window)) / step
                )
                overs.append(
                    100.0
                    * (max(p[1] for p in peak_window) - fast_level)
                    / step
                )
                tail_span = (slow.end - slow.start) * SETTLED_FRACTION
                drifts.append(100.0 * abs(slope(slow_tail)) * tail_span / step)
            result = PassResult(advance, k_idx)
            result.step = statistics.fmean(steps)
            result.undershoot_pct = statistics.fmean(unders)
            result.overshoot_pct = statistics.fmean(overs)
            result.drift_pct = statistics.fmean(drifts)
            result.spread = (
                statistics.pstdev(unders) if len(unders) > 1 else 0.0
            )
            result.edges = len(unders)
            results.append(result)
        results.sort(key=lambda r: r.advance)
        return results

    def _write_csv(self, gcmd, plan, samples, t_start):
        log_file = self.printer.get_start_args().get("log_file")
        directory = os.path.dirname(log_file) if log_file else "/tmp"
        path = os.path.join(
            directory, "pa_calibrate_%s.csv" % (time.strftime("%Y%m%d_%H%M%S"),)
        )
        edges = [(m.end, m.label, m.advance) for m in plan]
        try:
            with open(path, "w", newline="") as handle:
                writer = csv.writer(handle)
                writer.writerow(
                    ["advance", "time", "counts", "force_g", "segment"]
                )
                index = 0
                for sample in samples:
                    when = sample[SampleStructure.TIME.value]
                    while index < len(edges) and when > edges[index][0]:
                        index += 1
                    if index >= len(edges):
                        label, advance = "tail", edges[-1][2]
                    else:
                        _end, label, advance = edges[index]
                    writer.writerow(
                        [
                            "%.6f" % (advance,),
                            "%.6f" % (when - t_start,),
                            sample[SampleStructure.COUNTS.value],
                            "%.3f" % (sample[SampleStructure.FORCE_G.value],),
                            label,
                        ]
                    )
        except OSError as err:
            gcmd.respond_info("Could not write %s: %s" % (path, err))
            return
        gcmd.respond_info("Raw force trace written to %s" % (path,))

    def _report(self, gcmd, results, threshold):
        gcmd.respond_info(
            "advance  undershoot%  +/-   overshoot%  drift%     step  edges"
        )
        for res in results:
            gcmd.respond_info(
                "%7.4f  %10.2f %5.2f   %9.2f  %6.2f %8.0f %6d"
                % (
                    res.advance,
                    res.undershoot_pct,
                    res.spread,
                    res.overshoot_pct,
                    res.drift_pct,
                    res.step,
                    res.edges,
                )
            )
        unsettled = [r for r in results if r.drift_pct > MAX_TAIL_DRIFT * 100.0]
        if unsettled:
            gcmd.respond_info(
                "Verdict: INCONCLUSIVE. %i of %i passes never settled (tail "
                "drifts more than %.0f%% of the step over its own length), so "
                "their levels mean nothing. Raise PRIME and CYCLES."
                % (len(unsettled), len(results), MAX_TAIL_DRIFT * 100.0)
            )
            return
        # Below the knee the undershoot is flat but not at zero: the dip
        # metric carries a small per-run offset. Measure that flat part
        # and work relative to it, otherwise an absolute threshold moves
        # the answer whenever the offset does.
        flat = max(3, len(results) // 3)
        baseline = statistics.median(
            [r.undershoot_pct for r in results[:flat]]
        )
        level = baseline + threshold
        gcmd.respond_info(
            "Flat baseline over the lowest %i values: %.2f%%, so the knee is "
            "where undershoot passes %.2f%%." % (flat, baseline, level)
        )
        over = [r for r in results if r.undershoot_pct >= level]
        if not over:
            gcmd.respond_info(
                "Undershoot never rises %.1f%% above baseline. The value is "
                "above STOP=%.4f - extend the range."
                % (threshold, results[-1].advance)
            )
            return
        first = over[0]
        index = results.index(first)
        if index == 0:
            gcmd.respond_info(
                "Already %.2f%% above baseline at the lowest value tested. "
                "The value is below START=%.4f - extend the range downward."
                % (first.undershoot_pct - baseline, results[0].advance)
            )
            return
        # Interpolate across the bracketing pair so the answer is not
        # pinned to the sweep grid.
        low = results[index - 1]
        rise = first.undershoot_pct - low.undershoot_pct
        chosen = low.advance
        if rise > 0.0:
            chosen = low.advance + (level - low.undershoot_pct) * (
                first.advance - low.advance
            ) / rise
        gcmd.respond_info(
            "Knee is between %.4f (%.2f%%) and %.4f (%.2f%%), crossing at "
            "%.4f."
            % (
                low.advance,
                low.undershoot_pct,
                first.advance,
                first.undershoot_pct,
                chosen,
            )
        )
        gcmd.respond_info(
            "Set it with: SET_PRESSURE_ADVANCE ADVANCE=%.4f" % (chosen,)
        )


def load_config(config):
    return PressureAdvanceCalibrate(config)
