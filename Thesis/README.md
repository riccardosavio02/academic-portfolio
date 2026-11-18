# MSc thesis: impact-damaged thermal protection systems

MSc Aerospace Engineering (Aerogasdynamics), **Politecnico di Torino**, carried out at the **Université libre de Bruxelles (ULB)**
and the **von Karman Institute for Fluid Dynamics (VKI)** within the HITS project.
Supervisors: Prof. Thierry Magin (ULB/VKI), Prof. Domenic D'Ambrosio (Politecnico di Torino).

> *Thermomechanical behaviour and crack propagation of reusable Thermal Protection Systems subjected to high-velocity ballistic impact*

A piece of orbital debris hits a spacecraft at several km/s. What does the impact do to the shield, and how does the damaged
surface behave when it later faces a supersonic, high-enthalpy flow?

```
hypervelocity impact  ──►  damage in the shield  ──►  supersonic high-enthalpy flow over the damaged surface
  Abaqus/Explicit + VUMAT        cavity shape            Plasmatron tests  ◄──►  coupled flow–thermal CFD
        (ULB)                                                         (VKI)
```

| Folder | Content | Status |
|---|---|---|
| [ULB_Hypervelocity_Impact](./ULB_Hypervelocity_Impact) | User material (Fortran VUMAT), full model set-up to reproduce the simulations, results against LS-DYNA and experiments | Completed |
| [VKI_Plasmatron](./VKI_Plasmatron) | Aerothermodynamics of the damaged surface | In progress |
