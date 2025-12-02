# EFST-Master

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
![MATLAB](https://img.shields.io/badge/MATLAB-R2020b%2B-orange)
![Psychtoolbox](https://img.shields.io/badge/Psychtoolbox-3.0.18%2B-blue)
[![Built on respirationCA](https://img.shields.io/badge/built%20on-respirationCA-purple)](#foundation--attribution)

### Emotional & Color Stroop Toolbox for MATLAB / Psychtoolbox


---

## Overview

**EFST-Master** is a modular MATLAB/Psychtoolbox toolbox for running:

* **Color Stroop (3 colors)**
* **Emotional Face Stroop – 2 emotions (Happy / Sad)**
* **Emotional Face Stroop – 3 emotions (Happy / Neutral / Sad)**

It supports frame-accurate stimulus presentation, balanced trial sequencing, response logging, EEG triggers, and optional Pupil Labs eyetracking.

This toolbox is based on a **setup → sequence → run → post-processing → save** pipeline, designed for research-grade timing accuracy.

---

## Foundation & Attribution

This toolbox **builds upon, extends, and is structurally inspired by** the outstanding open-source framework:

### respirationCA Toolbox

**Martin Grund***, Esra Al, Marc Pabst, Alice Dabbagh, Tilman Stephani, Till Nierhaus,
Michael Gaebler, Arno Villringer

Max Planck Institute for Human Cognitive and Brain Sciences
Charité — Universitätsmedizin Berlin
Freie Universität Berlin

*Corresponding author: **[mgrund@cbs.mpg.de](mailto:mgrund@cbs.mpg.de)**

The code architecture, timing routines, sequence structure, and saving logic in EFST-Master originate from the foundational design of respirationCA, adapted and expanded for:

* emotional face paradigms
* 3-state congruency (C / I / IS)
* multi-emotion response options
* FACES database integration
* custom Stroop-style word–face pairing
* multi-task launch interface

This project would not exist without the respirationCA foundation.

---

## Emotional Face Stimuli — FACES Database

The Emotional Face Stroop tasks are compatible with the **FACES** image set:

**FACES Database**
Developed 2005–2007 by
Natalie C. Ebner, Michaela Riediger, Ulman Lindenberger
Center for Lifespan Psychology
Max Planck Institute for Human Development, Berlin

Website: [https://faces.mpdl.mpg.de/imeji/](https://faces.mpdl.mpg.de/imeji/)

Recommended filename structure:

```
999_y_f_h.jpg   % ID, age (young), female, happy
998_m_m_s.jpg   % ID, age (middle), male, sad
997_o_f_n.jpg   % ID, age (old),    female, neutral
```

The toolbox automatically parses: **age, gender, emotion**.

---

## Features

### Supported Tasks

* ✓ Color Stroop (3 ink colors)
* ✓ Emotional Face Stroop — **2 emotions**
* ✓ Emotional Face Stroop — **3 emotions**

### Timing & Display

* Frame-perfect timing using Psychtoolbox
* Photodiode patch for EEG synchronization
* Scaled face stimuli with centered overlay words
* ISI control (min, max, mean)
* fix → cue → stim → response pipeline

### Trial Sequencing

* Automatic parsing of FACES filenames
* Balanced:

  * age groups (y/m/o)
  * gender (f/m)
  * emotions (h/n/s)
  * congruency states (C / I / IS)
* Run-length avoidance for:

  * gender
  * age
  * emotion
  * congruency transitions
* No stimulus repetition across the entire experiment

### Response System

* 2-response version (Happy / Sad)
* 3-response version (Happy / Neutral / Sad)
* Customizable keyboard mappings
* Button box compatibility

### Eyetracking (optional)

* Pupil Labs (ZMQ) integration
* Sends annotations:

  * TrialStart
  * StimOnset
  * Response

### EEG Support

* Parallel-port triggers for:

  * trial start (30)
  * cue onset (11)
  * stimulus onset (12)
  * response screen onset (13)
  * response codes (40–42)

### Output

Each block produces a MATLAB `.mat` file containing:

* trial structure & metadata
* onsets (fixation, cue, stim, response)
* reaction times
* responses
* congruency
* age/gender/emotion of face
* ISIs
* EEG triggers (if enabled)

---

## Requirements

* MATLAB R2020b or newer
* Psychtoolbox 3.0.18+
* Parallel I/O driver (if EEG used)
* ZMQ toolbox (for Pupil Labs)
* FACES stimuli (or similar images)

---

## How to Run

1. Install Psychtoolbox.
2. Place your face images into `stimuli/`.
3. Open MATLAB in the EFST-Master directory.
4. Run:

```
run_master
```

5. Select:

   * task type
   * number of blocks
   * trials per block
   * EEG / eyetracking options
   * display screen

6. Follow on-screen instructions.

---

## Pipeline

```
setup → seq → run → intervals → save
```

**setup_xxx.m**
Initializes parameters (timing, keys, fonts, stimulus size).

**seq_xxx.m**
Builds balanced trial lists (age × gender × emotion × congruency).

**run_exp_xxx.m**
Runs the experiment (fix, cue, face, word, response).

**intervals.m**
Computes timing differences for data analysis.

**save_exp.m**
Saves behavioral + timing data to participant folder.

---

## Output Files

Per block:

```
faceStroop_3em_block_01_<ID>.mat
```

Global settings + full sequence:

```
faceStroop_3em_settings_seq.mat
```

---

## Citation Notice

If you publish work using EFST-Master, please cite:

### respirationCA Toolbox (foundational codebase)

**Grund, M., Al, E., Pabst, M., Dabbagh, A., Stephani, T., Nierhaus, T., Gaebler, M., & Villringer, A.**
Max Planck Institute for Human Cognitive and Brain Sciences
Contact: [mgrund@cbs.mpg.de](mailto:mgrund@cbs.mpg.de)

### FACES Database (stimulus set)

**Ebner, N. C., Riediger, M., & Lindenberger, U.**
Center for Lifespan Psychology
Max Planck Institute for Human Development, Berlin
[https://faces.mpdl.mpg.de/imeji/](https://faces.mpdl.mpg.de/imeji/)

### EFST-Master Toolbox

(You may cite this as a custom methodological tool in your dissertation/manuscript.)

---

## Contact

For questions, adjustments, or collaboration:

**Alparslan Önder**
(alparslanonder@med.uni-goettingen.de)

---


