# Load cell extrusion force bench test
#
# Answers one question: does this load cell see melt back-pressure when
# extruding into free air, and is the signal above the sensor noise floor?
# If it is, load-cell based pressure advance / max flow calibration is
# feasible on this hardware. If it is not, the load cell only sees nozzle
# contact force and the idea does not port.
#
# This file may be distributed under the terms of the GNU GPLv3 license.
from __future__ import annotations

import statistics

from klippy.extras.load_cell.load_cell import (
    LoadCell,
    SampleStructure,
    select_column,
)

# Fraction dropped from each end of a segment before averaging, to skip
# the accel/decel ramp and the pressure transient.
TRIM_FRACTION = 0.25
MIN_SEGMENT_SAMPLES = 8
# Baseline must relax to within this fraction of the measured step.
MAX_DRIFT_FRACTION = 0.1
# A segment whose own noise is this big a fraction of the step is not
# a steady state measurement.
MAX_SEGMENT_NOISE_FRACTION = 0.25


class Segment:
    def __init__(self, label, start, end, flow):
        self.label = label
        self.start = start
        self.end = end
        self.flow = flow
        self.mean = 0.0
        self.stddev = 0.0
        self.count = 0

    def measure(self, samples, column):
        span = self.end - self.start
        lo = self.start + span * TRIM_FRACTION
        hi = self.end - span * TRIM_FRACTION
        rows = [s for s in samples if lo <= s[SampleStructure.TIME.value] <= hi]
        self.count = len(rows)
        if self.count < MIN_SEGMENT_SAMPLES:
            return False
        values = select_column(rows, column)
        self.mean = statistics.fmean(values)
        self.stddev = statistics.pstdev(values)
        return True


class LoadCellFlowTest:
    def __init__(self, config):
        self.printer = config.get_printer()
        gcode = self.printer.lookup_object("gcode")
        gcode.register_command(
            "LOAD_CELL_FLOW_TEST",
            self.cmd_LOAD_CELL_FLOW_TEST,
            desc=self.cmd_LOAD_CELL_FLOW_TEST_help,
        )

    @staticmethod
    def _as_load_cell(obj) -> LoadCell | None:
        # [load_cell] registers a LoadCell directly. [load_cell_probe]
        # registers a LoadCellPrinterProbe that owns one privately, and
        # exposes no accessor for it.
        if isinstance(obj, LoadCell):
            return obj
        for attr in ("get_load_cell", "_load_cell"):
            inner = getattr(obj, attr, None)
            if callable(inner):
                inner = inner()
            if isinstance(inner, LoadCell):
                return inner
        return None

    def _find_load_cell(self, gcmd) -> tuple[str, LoadCell]:
        wanted = gcmd.get("LOAD_CELL", None)
        found: dict[str, LoadCell] = {}
        for obj_name, obj in self.printer.lookup_objects():
            load_cell = self._as_load_cell(obj)
            if load_cell is not None:
                found[obj_name] = load_cell
        if not found:
            raise gcmd.error(
                "No load cell found. Configure [load_cell] or "
                "[load_cell_probe]."
            )
        if wanted is None:
            if len(found) > 1:
                raise gcmd.error(
                    "Multiple load cells found, pass LOAD_CELL=<name>: %s"
                    % (", ".join(sorted(found)))
                )
            return next(iter(found.items()))
        for obj_name, load_cell in found.items():
            if wanted in (obj_name, obj_name.split()[-1]):
                return obj_name, load_cell
        raise gcmd.error(
            "Unknown load cell '%s', known: %s"
            % (wanted, ", ".join(sorted(found)))
        )

    cmd_LOAD_CELL_FLOW_TEST_help = (
        "Extrude into free air and report load cell force vs flow rate"
    )

    def cmd_LOAD_CELL_FLOW_TEST(self, gcmd):
        lc_name, load_cell = self._find_load_cell(gcmd)
        toolhead = self.printer.lookup_object("toolhead")
        gcode = self.printer.lookup_object("gcode")
        extruder = toolhead.get_extruder()
        if not extruder.get_name():
            raise gcmd.error("No extruder configured")
        flow_low = gcmd.get_float("FLOW_LOW", 2.0, above=0.0)
        flow_high = gcmd.get_float("FLOW_HIGH", 12.0, above=0.0)
        duration = gcmd.get_float("DURATION", 2.0, minval=0.5, maxval=10.0)
        if flow_high <= flow_low:
            raise gcmd.error("FLOW_HIGH must be greater than FLOW_LOW")
        # mm^3/s -> filament mm/min
        area = extruder.filament_area
        feed_low = (flow_low / area) * 60.0
        feed_high = (flow_high / area) * 60.0
        len_low = (flow_low / area) * duration
        len_high = (flow_high / area) * duration
        limit = extruder.max_e_dist
        for name, dist in (("FLOW_LOW", len_low), ("FLOW_HIGH", len_high)):
            if dist > limit:
                raise gcmd.error(
                    "%s needs a %.1fmm extrude-only move, over the %.1fmm "
                    "max_extrude_only_distance. Lower DURATION or the flow "
                    "rate." % (name, dist, limit)
                )
        total_mm = 2.0 * len_low + len_high
        gcmd.respond_info(
            "Flow test on '%s': idle/%.1f/%.1f/%.1f mm^3/s, %.1fs each, "
            "%.1fmm of filament. Extruding into free air - clear the nozzle."
            % (lc_name, flow_low, flow_high, flow_low, duration, total_mm)
        )

        # Start sampling before any of the test moves are queued.
        collector = load_cell.get_collector()
        t0 = toolhead.get_last_move_time()
        collector.start_collecting(min_time=t0)
        segments = []

        def queue(label, flow, script=None):
            start = toolhead.get_last_move_time()
            if script is None:
                toolhead.dwell(duration)
            else:
                gcode.run_script_from_command(script)
            end = toolhead.get_last_move_time()
            segments.append(Segment(label, start, end, flow))

        gcode.run_script_from_command("SAVE_GCODE_STATE NAME=_lc_flow_test")
        try:
            gcode.run_script_from_command("M83")
            queue("idle", 0.0)
            queue("low", flow_low, "G1 E%.4f F%.1f" % (len_low, feed_low))
            queue("high", flow_high, "G1 E%.4f F%.1f" % (len_high, feed_high))
            queue("low2", flow_low, "G1 E%.4f F%.1f" % (len_low, feed_low))
            queue("idle2", 0.0)
        finally:
            gcode.run_script_from_command(
                "RESTORE_GCODE_STATE NAME=_lc_flow_test"
            )
        samples, errors = collector.collect_until(segments[-1].end)
        self._report(gcmd, load_cell, segments, samples, errors)

    def _report(self, gcmd, load_cell, segments, samples, errors):
        if errors:
            gcmd.respond_info(
                "Sensor reported %i errors, %i overflows - results suspect"
                % (errors[0], errors[1])
            )
        if not samples:
            raise gcmd.error("No samples returned from sensor!")
        load_cell.validate_samples(samples, errors)
        column = SampleStructure.COUNTS.value
        for seg in segments:
            if not seg.measure(samples, column):
                raise gcmd.error(
                    "Only %i samples in the '%s' segment, need %i. Raise "
                    "DURATION or use a faster sensor."
                    % (seg.count, seg.label, MIN_SEGMENT_SAMPLES)
                )
        by_label = {seg.label: seg for seg in segments}
        idle = by_label["idle"]
        cpg = load_cell.get_counts_per_gram()

        def as_force(counts):
            if cpg is None:
                return ""
            return " (%+.2fg)" % (load_cell.invert * counts / cpg)

        gcmd.respond_info("flow      n    mean delta (counts)  noise (counts)")
        for seg in segments:
            delta = seg.mean - idle.mean
            gcmd.respond_info(
                "%-6s %5i  %+14.1f%s  %8.1f"
                % (seg.label, seg.count, delta, as_force(delta), seg.stddev)
            )
        noise = idle.stddev
        drift = abs(by_label["idle2"].mean - idle.mean)
        high_delta = abs(by_label["high"].mean - idle.mean)
        step_delta = abs(by_label["high"].mean - by_label["low"].mean)
        if noise <= 0.0:
            raise gcmd.error("Idle noise is zero - sensor is not streaming")
        gcmd.respond_info(
            "Baseline noise (1 sigma): %.1f counts%s"
            % (noise, as_force(abs(noise)))
        )
        gcmd.respond_info(
            "Baseline drift idle->idle2: %.1f counts (%.1f sigma)"
            % (drift, drift / noise)
        )
        gcmd.respond_info(
            "High flow vs idle:  %.1f sigma" % (high_delta / noise,)
        )
        gcmd.respond_info(
            "High flow vs low flow: %.1f sigma" % (step_delta / noise,)
        )
        snr = step_delta / noise
        # Sanity gates. Melt back-pressure relaxes when the extruder stops
        # and is steady within a segment. Anything else means the cell is
        # picking up something mechanical - usually extrudate fouling the
        # nozzle or the bed - and the separation figure is meaningless.
        problems = []
        if drift > MAX_DRIFT_FRACTION * step_delta:
            problems.append(
                "baseline did not return (%.0f%% of the measured step is "
                "still there after extrusion stopped) - the extrudate is "
                "probably touching the bed or the nozzle; park Z higher"
                % (100.0 * drift / step_delta,)
            )
        noise_limit = MAX_SEGMENT_NOISE_FRACTION * step_delta
        for seg in segments:
            if seg.flow and seg.stddev > noise_limit:
                problems.append(
                    "'%s' never reached steady state (noise %.0f%% of the "
                    "step)" % (seg.label, 100.0 * seg.stddev / step_delta)
                )
        low_delta = abs(by_label["low"].mean - idle.mean)
        if low_delta > high_delta:
            problems.append(
                "low flow loaded the cell more than high flow - force is not "
                "tracking flow rate"
            )
        if problems:
            gcmd.respond_info(
                "Verdict: INCONCLUSIVE (%.1f sigma separation, but the data "
                "is not trustworthy):" % (snr,)
            )
            for problem in problems:
                gcmd.respond_info("  - " + problem)
            return
        if snr >= 20.0:
            verdict = (
                "STRONG. This load cell sees melt back-pressure. Load-cell "
                "pressure advance tuning is worth building here."
            )
        elif snr >= 5.0:
            verdict = (
                "MARGINAL. There is a real signal but it is close to the "
                "noise. Try a larger FLOW_HIGH, a longer DURATION, or check "
                "for a hotend fan vibrating the toolhead."
            )
        else:
            verdict = (
                "WEAK. This load cell probably only sees nozzle contact "
                "force, not extruder reaction force. Free-air tuning will "
                "not work on this mounting."
            )
        gcmd.respond_info("Verdict: %.1f sigma separation. %s" % (snr, verdict))


def load_config(config):
    return LoadCellFlowTest(config)
