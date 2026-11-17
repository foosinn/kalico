# Measure what the load cell sees from toolhead acceleration alone
#
# The pressure advance sweep assumes acceleration is zero inside every
# measured window, because the cell reads toolhead inertia directly. An
# adaptive pressure advance table needs PA as a function of print
# acceleration, which breaks that assumption. This moves the toolhead at
# a range of accelerations with no extrusion at all and reports how much
# force each one produces, so the inertial contribution can be compared
# against the flow step the sweep measures.
#
# Nothing is extruded, so this costs no filament.
#
# This file may be distributed under the terms of the GNU GPLv3 license.
from __future__ import annotations

import csv
import os
import statistics
import time

from klippy.extras.load_cell.load_cell import SampleStructure
from klippy.extras.load_cell_flow_test import find_load_cell


class Leg:
    def __init__(self, accel, index):
        self.accel = accel
        self.index = index
        self.start = 0.0
        self.end = 0.0
        self.dwell = False


class LoadCellAccelTest:
    def __init__(self, config):
        self.printer = config.get_printer()
        gcode = self.printer.lookup_object("gcode")
        gcode.register_command(
            "LOAD_CELL_ACCEL_TEST",
            self.cmd_LOAD_CELL_ACCEL_TEST,
            desc=self.cmd_LOAD_CELL_ACCEL_TEST_help,
        )

    cmd_LOAD_CELL_ACCEL_TEST_help = (
        "Move without extruding and report load cell force per acceleration"
    )

    def cmd_LOAD_CELL_ACCEL_TEST(self, gcmd):
        lc_name, load_cell = find_load_cell(self.printer, gcmd)
        toolhead = self.printer.lookup_object("toolhead")
        gcode = self.printer.lookup_object("gcode")
        reactor = self.printer.get_reactor()
        status = toolhead.get_status(reactor.monotonic())
        if not all(a in status.get("homed_axes", "") for a in "xyz"):
            raise gcmd.error("Home all axes first")

        raw = gcmd.get("ACCELS", "500,1000,2000,4000")
        try:
            accels = [float(v) for v in raw.split(",") if v.strip()]
        except ValueError:
            raise gcmd.error("ACCELS must be a comma separated list")
        if not accels:
            raise gcmd.error("ACCELS is empty")
        limit = status.get("max_accel")
        over = [a for a in accels if limit and a > limit]
        if over:
            raise gcmd.error(
                "ACCELS %s exceed the machine max_accel of %.0f"
                % (", ".join("%.0f" % a for a in over), limit)
            )
        speed = gcmd.get_float("SPEED", 100.0, above=0.0)
        distance = gcmd.get_float("DISTANCE", 60.0, above=0.0)
        repeats = gcmd.get_int("REPEATS", 4, minval=1, maxval=50)
        settle = gcmd.get_float("SETTLE", 0.5, minval=0.2, maxval=5.0)
        axis = gcmd.get("AXIS", "X").strip().upper()
        if axis not in ("X", "Y"):
            raise gcmd.error("AXIS must be X or Y")
        # The slowest acceleration needs room for a full ramp each way,
        # otherwise it never reaches the speed and the test understates
        # what that acceleration does.
        ramp = speed * speed / (2.0 * min(accels))
        if 2.0 * ramp > distance:
            raise gcmd.error(
                "DISTANCE %.1fmm is too short: %.0f mm/s at %.0f mm/s^2 "
                "needs %.1fmm of ramp each way. Raise DISTANCE or SPEED "
                "down." % (distance, speed, min(accels), ramp)
            )

        gcode_move = self.printer.lookup_object("gcode_move")
        gc_pos = gcode_move.get_status()["gcode_position"]
        th_pos = toolhead.get_position()
        x0 = gcmd.get_float("X", gc_pos[0])
        y0 = gcmd.get_float("Y", gc_pos[1])
        z0 = gcmd.get_float("Z", gc_pos[2])
        index = 0 if axis == "X" else 1
        offset = th_pos[index] - gc_pos[index]
        base = x0 if axis == "X" else y0
        axis_max = status.get("axis_maximum")
        if axis_max is not None and base + distance + offset > axis_max[index]:
            raise gcmd.error(
                "Test needs %s up to %.1f, past the %.1f limit. Lower "
                "DISTANCE or pass a smaller %s."
                % (axis, base + distance + offset, axis_max[index], axis)
            )

        old = (
            status["max_velocity"],
            status["max_accel"],
            status["square_corner_velocity"],
            status["minimum_cruise_ratio"],
        )
        gcmd.respond_info(
            "Acceleration test on '%s': %s mm/s^2, %.0f mm/s over %.0fmm in "
            "%s, %i repeats each. No extrusion."
            % (
                lc_name,
                ", ".join("%.0f" % a for a in accels),
                speed,
                distance,
                axis,
                repeats,
            )
        )

        gcode.run_script_from_command("SAVE_GCODE_STATE NAME=_lc_accel_test")
        try:
            gcode.run_script_from_command("G90")
            samples, idle_end, t_start, legs = self._run(
                gcmd, load_cell, toolhead, gcode, accels, speed, distance,
                repeats, settle, axis, base, x0, y0, z0
            )
        finally:
            gcode.run_script_from_command(
                "SET_VELOCITY_LIMIT VELOCITY=%.6f ACCEL=%.6f "
                "SQUARE_CORNER_VELOCITY=%.6f MINIMUM_CRUISE_RATIO=%.6f"
                % old
            )
            gcode.run_script_from_command(
                "RESTORE_GCODE_STATE NAME=_lc_accel_test MOVE=0"
            )
        self._write_csv(gcmd, legs, samples, t_start)
        self._report(gcmd, load_cell, legs, samples, idle_end, t_start)

    def _run(self, gcmd, load_cell, toolhead, gcode, accels, speed, distance,
             repeats, settle, axis, base, x0, y0, z0):
        feed = speed * 60.0
        gcode.run_script_from_command(
            "G1 X%.3f Y%.3f Z%.3f F%.0f" % (x0, y0, z0, 3000.0)
        )
        collector = load_cell.get_collector()
        t_start = toolhead.get_last_move_time()
        collector.start_collecting(min_time=t_start)
        toolhead.dwell(settle)
        idle_end = toolhead.get_last_move_time()

        legs = []

        def mark(leg):
            def callback(print_time):
                leg.end = print_time

            toolhead.register_lookahead_callback(callback)

        for accel in accels:
            # Full trapezoid, so the commanded acceleration is actually
            # reached rather than being clipped by the cruise ratio.
            gcode.run_script_from_command(
                "SET_VELOCITY_LIMIT ACCEL=%.6f MINIMUM_CRUISE_RATIO=0"
                % (accel,)
            )
            for repeat in range(repeats):
                for target in (base + distance, base):
                    leg = Leg(accel, len(legs))
                    gcode.run_script_from_command(
                        "G1 %s%.4f F%.0f" % (axis, target, feed)
                    )
                    mark(leg)
                    legs.append(leg)
            # Let the structure ring down before the next acceleration
            toolhead.dwell(settle)
            legs.append(Leg(accel, len(legs)))
            legs[-1].end = toolhead.get_last_move_time()
            legs[-1].dwell = True

        t_end = toolhead.get_last_move_time()
        samples, errors = collector.collect_until(t_end)
        if not samples:
            raise gcmd.error("No samples returned from sensor!")
        load_cell.validate_samples(samples, errors)
        previous = idle_end
        for leg in legs:
            if not leg.end:
                raise gcmd.error("A lookahead callback never fired")
            leg.start = previous
            previous = leg.end
        return samples, idle_end, t_start, legs

    def _write_csv(self, gcmd, legs, samples, t_start):
        log_file = self.printer.get_start_args().get("log_file")
        directory = os.path.dirname(log_file) if log_file else "/tmp"
        path = os.path.join(
            directory,
            "accel_test_%s.csv" % (time.strftime("%Y%m%d_%H%M%S"),),
        )
        try:
            with open(path, "w", newline="") as handle:
                writer = csv.writer(handle)
                writer.writerow(
                    ["accel", "time", "counts", "force_g", "phase"]
                )
                index = 0
                for sample in samples:
                    when = sample[SampleStructure.TIME.value]
                    while index < len(legs) and when > legs[index].end:
                        index += 1
                    if index >= len(legs):
                        accel, phase = 0.0, "tail"
                    else:
                        leg = legs[index]
                        accel = leg.accel
                        phase = "dwell" if leg.dwell \
                            else "move"
                    writer.writerow(
                        [
                            "%.0f" % (accel,),
                            "%.6f" % (when - t_start,),
                            sample[SampleStructure.COUNTS.value],
                            "%.3f" % (sample[SampleStructure.FORCE_G.value],),
                            phase,
                        ]
                    )
        except OSError as err:
            gcmd.respond_info("Could not write %s: %s" % (path, err))
            return
        gcmd.respond_info("Raw force trace written to %s" % (path,))

    def _report(self, gcmd, load_cell, legs, samples, idle_end, t_start):
        column = SampleStructure.COUNTS.value
        idle = [
            s[column]
            for s in samples
            if t_start <= s[SampleStructure.TIME.value] <= idle_end
        ]
        if len(idle) < 2:
            raise gcmd.error("Too few idle samples for a noise estimate")
        base = statistics.fmean(idle)
        noise = statistics.pstdev(idle)
        cpg = load_cell.get_counts_per_gram()
        by_accel = {}
        for leg in legs:
            if leg.dwell:
                continue
            window = [
                s[column]
                for s in samples
                if leg.start <= s[SampleStructure.TIME.value] <= leg.end
            ]
            if not window:
                continue
            peak = max(abs(v - base) for v in window)
            by_accel.setdefault(leg.accel, []).append(peak)
        gcmd.respond_info(
            "Idle noise 1 sigma: %.0f counts%s"
            % (noise, "" if cpg is None else " (%.1fg)" % (noise / cpg,))
        )
        gcmd.respond_info("accel     peak deviation      sigma   legs")
        for accel in sorted(by_accel):
            peaks = by_accel[accel]
            peak = statistics.fmean(peaks)
            grams = "" if cpg is None else " (%7.1fg)" % (peak / cpg,)
            gcmd.respond_info(
                "%6.0f  %10.0f counts%s  %7.1f %6d"
                % (accel, peak, grams, peak / noise if noise else 0.0,
                   len(peaks))
            )
        gcmd.respond_info(
            "Compare these against the flow step CALIBRATE_PRESSURE_ADVANCE "
            "measures. If they are a small fraction of it, an acceleration "
            "axis is workable; if not, the inertial signal has to be "
            "subtracted with a matching dry run."
        )


def load_config(config):
    return LoadCellAccelTest(config)
