# Load cell pressure advance calibration

This document describes calibrating [pressure
advance](Pressure_Advance.md) from load cell force data, without
printing or judging a test object.

The idea comes from [CNC
Kitchen](https://www.cnckitchen.com/blog/i-automated-pressure-advance-using-a-bed-leveling-sensor),
who used the load cell in a Prusa print head for the same purpose. It
requires a [load cell](Load_Cell.md) that senses extruder reaction
force, which is the case when the cell sits under the hotend or
extruder assembly. A cell mounted only to read nozzle contact force in
Z will not work. `LOAD_CELL_FLOW_TEST` below answers that question for
a given machine before anything else is attempted.

Warning: these commands extrude several grams of filament into free
air. Park the nozzle high above the bed, clear of anything it could
touch.

## How it works

While extruding, the force on the load cell tracks the pressure in the
melt zone. Step the flow rate from slow to fast and back, and the force
follows with a lag. With too little pressure advance the force rounds
off at each transition. With too much, it overshoots on the way up and
dips below the settled level on the way down. That dip, the
*undershoot*, is the signal: it is absent below the correct value and
grows steadily above it. The calibration sweeps pressure advance and
reports where the undershoot starts.

Three properties of Kalico's motion system make this practical:

- Force samples are timestamped on the same clock as the move queue, so
  each flow transition is located exactly rather than estimated.
- `toolhead.register_lookahead_callback()` returns the print time at a
  move boundary without flushing the queue, so the segments blend and
  the flow step stays sharp.
- Pressure advance can be changed between moves without stopping
  motion, as long as `pressure_advance_smooth_time` does not change.

## Checking the load cell first

```
LOAD_CELL_FLOW_TEST FLOW_LOW=2 FLOW_HIGH=24 DURATION=2
```

Extrudes at idle, low, high and low flow in turn while stationary, then
reports the force in each segment. Add `[load_cell_flow_test]` to the
config to enable it.

What matters is the verdict line. A strong separation means the cell
sees melt back-pressure and the calibration below is worth running.
`FLAT` with filament loaded means it does not, and the method does not
apply to that machine.

The test refuses to give a verdict when the data is untrustworthy: if
the baseline does not return after extrusion stops, or a segment never
reaches steady state, something mechanical is loading the cell.
Extrudate touching the bed is the usual cause, so park Z high.

## Calibrating pressure advance

```ini
[pressure_advance_calibrate]
```

```
G1 Z100 F600
CALIBRATE_PRESSURE_ADVANCE START=0.03 STOP=0.06 STEP=0.0025
```

The nozzle must be hot, filament loaded and all axes homed. A sweep of
13 values takes about five minutes and uses roughly 4g of filament.

Start with a coarse sweep (`START=0 STOP=0.08 STEP=0.01`) to find the
neighbourhood, then refine around the knee.

### Output

```
advance  undershoot%  +/-   overshoot%  drift%     step  edges
 0.0400        1.05  0.08        3.29    3.66   160787      2
 0.0475        1.17  0.37        3.87    6.20   164535      2
 0.0500        6.05  0.06        2.87    4.78   165421      2
 ...
Flat baseline over the lowest 4 values: -0.13%, so the knee is where
undershoot passes 2.87%.
Knee is between 0.0475 (1.17%) and 0.0500 (6.05%), crossing at 0.0484.
```

- `undershoot%` is the metric. It is flat below the knee and rises
  steadily above it.
- `+/-` is the spread between the measured edges. Large values mean the
  result is not repeatable.
- `overshoot%` rises with pressure advance as well, since it inflates
  both transitions. It is a cross-check, not a selection criterion.
- `drift%` is how far a settled window moves over its own length. Above
  10% the pass never settled and no verdict is given.
- `step` is the force difference between the flow rates, in ADC counts.

The raw force trace is written to a CSV beside `klippy.log` and can be
plotted with `scripts/graph_pa_calibrate.py`. Read it whenever a result
looks wrong; the aggregate numbers hide the shape.

### Useful parameters

- `FLOW_LOW`, `FLOW_HIGH`\
  The two flow rates, in mm³/s. Their difference sets the size of the
  signal, because pressure advance injects `pa` times the change in
  extruder velocity. The absolute flow rate does not matter nearly as
  much as the difference.

- `PRIME`\
  Seconds of cycling before measuring, default 30. Melt pressure needs
  far longer to equilibrate than one flow step, and priming matters
  more than any other parameter. See the note below.

- `CYCLES`\
  Cycles per value, default 5. The cycles after each change of pressure
  advance are discarded and the rest are averaged.

- `DURATION`\
  Seconds per segment, default 2. Long enough that the slow segment
  settles before the end, which is what the `drift%` column checks.

- `ALTERNATE`\
  `1` reverses travel direction per value and keeps the purge in a
  band. `0` lays one continuous line, which needs more bed width.

- `SEED`\
  The sweep runs values in a shuffled but reproducible order, because
  pressure advance rising with time is indistinguishable from drift.

## Notes on the measurement

**Pressure advance requires XY motion.** Kalico applies no pressure
advance to extrude-only moves, so every segment travels. Travel speed
is constant throughout and the flow step comes from changing extrusion
width, which keeps toolhead motion out of the force trace.

**Travel limits the flow step.** The lookahead caps junction speed
between the segments at `instantaneous_corner_velocity` divided by the
change in extrusion ratio. Travel speed cancels out, leaving a limit on
the flow difference alone. The command raises
`instantaneous_corner_velocity` for the duration of the sweep and
restores it afterwards, along with `max_extrude_cross_section`, which
slow travel would otherwise exceed.

**Melt pressure equilibrates slowly.** Starting and stopping extrusion
for each value measures the startup transient instead of pressure
advance, and that transient grows with pressure advance, which inverts
the metric. The sweep therefore extrudes continuously from the prime to
the end, and primes with the same cycle it will measure, since the melt
equilibrates to the mean flow rather than the low one.

**What the load cell responds to.** `LOAD_CELL_ACCEL_TEST` moves the
toolhead without extruding, at a range of accelerations. On the machine
this was developed against, a sixteenfold range of acceleration changed
the force by 7%, while position along the axis changed it by about
0.5g/mm. Acceleration is not the hazard; slow position-dependent drift
is, which is why passes are kept short and `drift%` is checked.

## Limitations

The knee is where over-compensation first becomes measurable, which is
not necessarily where prints look best. CNC Kitchen reported their
automatic values landing slightly above visually chosen ones, and the
same offset appeared here. Treat the result as a well-defined,
repeatable starting point and confirm it against a print.

Values repeat to about ±0.003 between runs. The result applies to the
filament, nozzle and temperature it was measured at, and on a tool
changer to the tool that was mounted.
