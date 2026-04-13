---
name: R coding style
description: User wants curly brackets on all loops and one command per line in R code
type: feedback
---

Always use curly brackets for loops (even one-liners) and put each command on its own line in R code.

**Why:** User prefers explicit, readable R code matching a MATLAB-like style — no compressed one-liners.
**How to apply:** Any `for`, `if`, `else` block gets `{ }` on separate lines, and never chain multiple statements with `;` on one line.
