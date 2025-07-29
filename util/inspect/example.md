<!--
SPDX-FileCopyrightText: 2006-2025, Knut Reinert & Freie Universität Berlin
SPDX-FileCopyrightText: 2016-2025, Knut Reinert & MPI für molekulare Genetik
SPDX-License-Identifier: CC-BY-4.0
-->

# Nomenclature

For an IBF, `User bins` and `Active bins` are the same.
In an HIBF,
  * `User bins`: Bins that actually contain the data of a user bin.
  * `Active bins`: Merged or deleted bins (see `seqan::hibf::bin_kind`).

# Format

* The data for each IBF is printed in a block.
* The children of an IBF are printed after the parent.
* Intendation levels are used to indicate the hierarchy of IBFs, e.g., the children of an IBF are indented by one level.
* A level of indentation is four spaces.
* The data for an IBF in the HIBF is indented by two spaces.
* The [HIBF with Format Annotations](#hibf-with-format-annotations) section shows an output with explanations/annotations.

# Miscellaneous

* The output is colored on the terminal.
* The output is not colored when redirected/piped.
* Test data may not be up-to-date if IBF/HIBF implementations have changed.
* Test data added: 2025-07-29

# IBF

**Example data:** https://ftp.seqan.de/hibf/util.test.ibf \
**Command:** `./inspect --input util.test.ibf`

```text
User bins: 128
Technical bins: 128
Empty bins: 0
Bin words: 2
Bin size: 1,024
Total size: 131,072
Hash functions: 3
Hash shift: 53
```

# Simple HIBF

**Example data:** https://ftp.seqan.de/hibf/util.test.simple.hibf \
**Command:** `./inspect --input util.test.simple.hibf --hibf`

```text
User bins: 2
Empty bins: 0
IBFs: 1
Total size: 4,112
ID: 0 Level: 0 Children: 0
  User bins: 2
  Active bins: 64
  Technical bins: 64
  Empty bins: 0
  Bin words: 1
  Bin size: 48
  Total size: 3,072
  Hash functions: 2
  Hash shift: 58
```

# HIBF

**Example data:** https://ftp.seqan.de/hibf/util.test.hibf \
**Command:** `./inspect --input util.test.hibf --hibf`

```text
User bins: 4,097
Empty bins: 0
IBFs: 66
Total size: 5,362,464
ID: 0  Level: 0 Children: 64
  User bins: 0
  Active bins: 64
  Technical bins: 64
  Empty bins: 0
  Bin words: 1
  Bin size: 2,642
  Total size: 169,088
  Hash functions: 2
  Hash shift: 52
    ID: 3  Level: 1 Children: 0
      User bins: 64
      Active bins: 64
      Technical bins: 64
      Empty bins: 0
      Bin words: 1
      Bin size: 1,245
      Total size: 79,680
      Hash functions: 2
      Hash shift: 53
    [...]
    ID: 1  Level: 1 Children: 1
      User bins: 63
      Active bins: 64
      Technical bins: 64
      Empty bins: 0
      Bin words: 1
      Bin size: 1,245
      Total size: 79,680
      Hash functions: 2
      Hash shift: 53
        ID: 2  Level: 2 Children: 0
          User bins: 2
          Active bins: 64
          Technical bins: 64
          Empty bins: 0
          Bin words: 1
          Bin size: 394
          Total size: 25,216
          Hash functions: 2
          Hash shift: 55
```

# HIBF with Format Annotations

```text
User bins: 4,097 ──────┐
Empty bins: 0          ├ HIBF Metadata
IBFs: 66               │
Total size: 5,362,464 ─┘
ID: 0  Level: 0 Children: 64 ─┬── Root IBF
  User bins: 0                │
  Active bins: 64             │
  Technical bins: 64          │
  Empty bins: 0               │
  Bin words: 1                ├ Root IBF Metadata
  Bin size: 2,642             │
  Total size: 169,088         │
  Hash functions: 2           │
  Hash shift: 52 ─────────────┘
    ID: 3  Level: 1 Children: 0 ──────┐
      User bins: 64                   │
      Active bins: 64                 │
      Technical bins: 64              │
      Empty bins: 0                   │
      Bin words: 1                    │
      Bin size: 1,245                 │
      Total size: 79,680              │
      Hash functions: 2               │
      Hash shift: 53                  │
    [...]                             ├ Direct Children if the Root IBF
    ID: 1  Level: 1 Children: 1       │
      User bins: 63                   │
      Active bins: 64                 │
      Technical bins: 64              │
      Empty bins: 0                   │
      Bin words: 1                    │
      Bin size: 1,245                 │
      Total size: 79,680              │
      Hash functions: 2               │
      Hash shift: 53 ─────────────────┘
        ID: 2  Level: 2 Children: 0 ─┐
          User bins: 2               │
          Active bins: 64            │
          Technical bins: 64         │
          Empty bins: 0              │
          Bin words: 1               ├ Child of IBF with ID=1
          Bin size: 394              │
          Total size: 25,216         │
          Hash functions: 2          │
          Hash shift: 55 ────────────┘
```
