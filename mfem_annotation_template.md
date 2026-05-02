# 📘 [Example Title]

`[duration]` `[difficulty: basic | intermediate | advanced]`

> **✓ Lesson Objectives**
>
> - [Objective 1 — what the reader will understand]
> - [Objective 2 — what the reader will be able to do]
> - [Objective 3 — optional]

> **ℹ Note**
>
> Please complete the [Finite Element Basics](https://mfem.org/tutorial/fem/) lesson on `ex1.cpp` before this one. [Add any other prerequisite examples here, e.g., "We also recommend viewing Example N before this one."]

---

## ☑ [Problem name]

[One-paragraph description of the PDE problem this example solves. State what the PDE models physically — heat conduction, electromagnetics, elasticity, eigenmodes, etc.]

[State the strong form as a numbered display equation:]

$$\text{[strong form of the PDE]} \tag{1}$$

[with boundary conditions:]

$$\text{[boundary conditions]}$$

[where ... explain each symbol that appears.]

### Weak form

[Describe the test-function multiplication and integration by parts, the same way the `ex1` tutorial does. Show the steps as numbered equations.]

Multiplying (1) by a test function $\varphi_i$ and integrating by parts:

$$\text{[after multiplication]} \tag{2}$$

$$\text{[after integration by parts]} \tag{3}$$

[Note any boundary terms that vanish, and why.]

Substituting the FE expansion $u_h = \sum_j c_j \varphi_j$ gives the matrix system

$$\text{[matrix form]} \tag{4}$$

where

$$A_{ij} = \text{[bilinear form entry]} \tag{5}$$

$$b_i = \text{[linear form entry]} \tag{6}$$

$$x_j = c_j \tag{7}$$

[Brief paragraph: properties of $A$ — symmetric? positive definite? saddle-point? — and what that implies for solvers.]

### [Optional subsection: anything specific to this example]

[E.g., "Why a generalized eigenvalue problem?", "Why an $H(\text{curl})$ space?", "Why DG?". One short subsection per non-obvious choice the example makes.]

> **ℹ Note**
>
> [Optional pointer to deeper theory or related MFEM docs. Keep it brief — one or two sentences.]

---

## ☑ Annotated Example [N]

MFEM's Example [N] implements the above formulation in the source file [`examples/ex[N]p.cpp`](https://github.com/mfem/mfem/blob/master/examples/ex[N]p.cpp). [State whether the example has a serial version, a parallel version, or both, and which you're annotating.]

[One-paragraph summary of what the example does end-to-end: e.g., "We compute the lowest `nev` eigenpairs on a mesh provided as input."]

Below we highlight selected portions of the example code and connect them with the description in the previous section. You can follow along by browsing [`ex[N]p.cpp`](https://github.com/mfem/mfem/blob/master/examples/ex[N]p.cpp) in your editor.

### [Section 1 — typically MPI/HYPRE init for a parallel example]

[One-line description of what this block does.] ([lines X–Y](https://github.com/mfem/mfem/blob/master/examples/ex[N]p.cpp#LX-LY)):

```cpp
[code excerpt]
```

[Prose explanation: what the code does, why it's here, anything subtle.]

### [Section 2 — command-line options]

[Short description.] ([lines X–Y](https://github.com/mfem/mfem/blob/master/examples/ex[N]p.cpp#LX-LY)):

```cpp
[code excerpt]
```

[Explain what each option controls. Highlight the ones specific to this example.]

### [Section 3 — mesh construction]

[Description.] ([lines X–Y](https://github.com/mfem/mfem/blob/master/examples/ex[N]p.cpp#LX-LY)):

```cpp
[code excerpt]
```

[Explain the serial→parallel→refine pattern, or whatever mesh handling is specific to this example.]

### [Section 4 — finite element space]

[Description.] ([lines X–Y](https://github.com/mfem/mfem/blob/master/examples/ex[N]p.cpp#LX-LY)):

```cpp
[code excerpt]
```

[Explain which FE space is being built (H1, H(curl), H(div), L2) and why it's the right space for this problem.]

### [Section 5 — boundary conditions]

[Description.] ([lines X–Y](https://github.com/mfem/mfem/blob/master/examples/ex[N]p.cpp#LX-LY)):

```cpp
[code excerpt]
```

[Explain how essential vs. natural BCs are handled. If there's anything tricky about the BC handling — e.g., elimination, weak imposition — explain it here.]

### [Section 6 — bilinear / linear forms]

[Description.] ([lines X–Y](https://github.com/mfem/mfem/blob/master/examples/ex[N]p.cpp#LX-LY)):

```cpp
[code excerpt]
```

[Connect each `Add...Integrator(...)` call back to the corresponding term in the weak form (equations 2–3 above). Mention which integrators are used and what they correspond to mathematically.]

[If there's anything subtle here — like a special trick the example uses to handle a singular operator, or a non-standard boundary integrator — call it out:]

**(i) [Subtlety name].** [Explanation.]

**(ii) [Another subtlety, if applicable].** [Explanation.]

### [Section 7 — preconditioner / solver setup]

[Description.] ([lines X–Y](https://github.com/mfem/mfem/blob/master/examples/ex[N]p.cpp#LX-LY)):

```cpp
[code excerpt]
```

[Explain the choice of solver and preconditioner. Why is this combination appropriate for this PDE? What's the expected scaling behavior?]

### [Section 8 — solve]

[Description.] ([lines X–Y](https://github.com/mfem/mfem/blob/master/examples/ex[N]p.cpp#LX-LY)):

```cpp
[code excerpt]
```

[Explain what `Solve()` / `Mult()` does, what the return value or output is, and how the solution is post-processed if needed.]

### [Section 9 — output / visualization]

[Description.] ([lines X–Y](https://github.com/mfem/mfem/blob/master/examples/ex[N]p.cpp#LX-LY)):

```cpp
[code excerpt]
```

[Explain how the solution is written to disk and/or sent to GLVis.]

[Add additional sections as needed — e.g., error computation against an exact solution, time-stepping loop, AMR loop. Use the same pattern: description → code excerpt → prose.]

---

## ☑ Sample runs

A few representative invocations (these match the comments at the top of `ex[N]p.cpp`):

```bash
mpirun -np 4 ex[N]p -m ../data/[mesh1].mesh
mpirun -np 4 ex[N]p -m ../data/[mesh2].mesh -o 2
mpirun -np 4 ex[N]p -m ../data/[mesh3].mesh -[other flag]
```

[Brief description of what each sample run is testing — different mesh types, different orders, different physics regimes.]

> **ℹ Try this!**
>
> [Concrete experiment 1 — usually about convergence or basic behavior. Include expected output or a property the reader can check.]
>
> ```
> [expected output snippet, if useful]
> ```

> **ℹ Try this!**
>
> [Concrete experiment 2 — usually exploring a parameter (polynomial order, mesh refinement, problem variant). State what theoretical prediction the reader should verify.]

> **⚠ Warning**
>
> [Optional warning about a common pitfall — e.g., "Don't forget the `./` before the executable on macOS", or "This example requires GLVis to be running on port 19916 to see visualizations".]

---

## ☑ [Optional: special section unique to this example]

[Use this slot for content that doesn't fit elsewhere. Examples:]

- **Why parallel-only?** — if the example has no serial version, explain why.
- **Convergence study** — if the example is well-suited to a discretization-error study.
- **Parallel scaling notes** — if the example demonstrates a performance feature.
- **Comparison to another example** — e.g., "Compared to ex1, this example differs in...".

[Delete this section if not needed.]

---

## ☑ Things to keep in mind

- **[Key takeaway 1 — usually the most important conceptual point].** [One- or two-sentence elaboration.]

- **[Key takeaway 2].** [Elaboration.]

- **[Key takeaway 3].** [Elaboration.]

- **[Add 1–2 more bullets if needed.]**

---

## ☑ Next Steps

- [Example M](https://mfem.org/examples/#exM): [One-line description of why a reader of this annotation might want to look at it next.]
- [Example K](https://mfem.org/examples/#exK): [Description.]
- [Optional: link to a miniapp or tutorial page for further reading.]

[Back to the MFEM tutorial page](https://mfem.org/tutorial/)

---

*Based on `ex[N]p.cpp` from MFEM master branch. Line numbers refer to the master version on GitHub at the time of writing and may shift slightly in other releases.*

*Annotated by [Your Name] and [Partner's Name], APMA 2560, [Date].*
