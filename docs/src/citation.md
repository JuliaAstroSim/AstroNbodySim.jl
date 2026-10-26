# Citing AstroNbodySim.jl

If you use `AstroNbodySim.jl` in published research, please cite the package
repository together with the exact version or commit you used (`]status` in
the Julia REPL prints both).

## Primary citation

`AstroNbodySim.jl` is the N-body engine of
[WaveDM.jl](https://github.com/JuliaAstroSim/WaveDM.jl), and was developed in
the course of the following work:

> Meng, R. & Dong, X. 2026, *WaveDM.jl: An Adaptable Simulation Framework for
> Dynamics of Baryonic and Wave Dark Matter on Galaxy Scales*,
> Research in Astronomy and Astrophysics.
> [Journal article](http://iopscience.iop.org/article/10.1088/1674-4527/ae9042)
> &middot; [arXiv:2606.25026](https://arxiv.org/abs/2606.25026)
> &middot; DOI [10.1088/1674-4527/ae9042](https://doi.org/10.1088/1674-4527/ae9042)

```bibtex
@article{meng2026wavedmjladaptablea,
  title = {{{WaveDM}}.Jl: {{An Adaptable Simulation Framework}} for {{Dynamics}} of {{Baryonic}} and {{Wave Dark Matter}} on {{Galaxy Scales}}},
  shorttitle = {{{WaveDM}}.Jl},
  author = {Meng, Runyu and Dong, Xiaobo},
  year = 2026,
  journal = {Research in Astronomy and Astrophysics},
  issn = {1674-4527},
  doi = {10.1088/1674-4527/ae9042},
  url = {http://iopscience.iop.org/article/10.1088/1674-4527/ae9042},
  eprint = {2606.25026},
  archiveprefix = {arXiv},
  primaryclass = {astro-ph.GA},
  urldate = {2026-07-26}
}
```

## Publications built on this code

The following papers were produced with this code base and should be cited in
addition to the entry above when they are the work your results refer to:

> Meng, R. & Dong, X. 2026, *Wave Cold Dark Matter Can Reproduce the Radial
> Acceleration Relation in X-ray Elliptical Galaxies*,
> Research in Astronomy and Astrophysics.
> [Journal article](http://iopscience.iop.org/article/10.1088/1674-4527/aeac08)
> &middot; [arXiv:2610.18969](https://arxiv.org/abs/2610.18969)
> &middot; DOI [10.1088/1674-4527/aeac08](https://doi.org/10.1088/1674-4527/aeac08)

```bibtex
@article{meng2026wavecold,
  title = {Wave {{Cold Dark Matter Can Reproduce}} the {{Radial Acceleration Relation}} in {{X}}-Ray {{Elliptical Galaxies}}},
  author = {Meng, Runyu and Dong, Xiaobo},
  year = 2026,
  journal = {Research in Astronomy and Astrophysics},
  issn = {1674-4527},
  doi = {10.1088/1674-4527/aeac08},
  url = {http://iopscience.iop.org/article/10.1088/1674-4527/aeac08},
  eprint = {2610.18969},
  archiveprefix = {arXiv},
  primaryclass = {astro-ph.GA},
  urldate = {2026-09-24}
}
```

## Methodological references

When citing a specific gravity solver or MOND model, also cite the original
methodological paper — for example:

- **DirectSum** — Aarseth, S. J. 1963, MNRAS, 126, 223.
- **Gadget-2 cubic-spline softening and leapfrog integration** — Springel, V., Yoshida, N. & White, S. D. M. 2001, New Astronomy, 6, 79; Springel, V. 2005, Nature, 435, 629.
- **Tree (octree) gravity** — Barnes, J. & Hut, P. 1986, Nature, 324, 446.
- **MOND (Milgrom 1983)** — Milgrom, M. 1983, ApJ, 270, 365.
- **QUMOND** — Milgrom, M. 2010, MNRAS, 403, 886.
- **Pseudo-Newtonian potential** — Paczyński, B. & Wiita, P. J. 1980, A&A, 88, 23.

