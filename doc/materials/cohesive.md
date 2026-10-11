# Cohesive Zone Models

Cohesive constitutive models describe the local relation between the
displacement jump across two adjoining surfaces and the tractions transferred
between them. The response generally includes an initial elastic or penalty
regime, followed by damage and softening as the interface opens or slides, and
may end in complete decohesion. Parameters such as the peak traction and
fracture energy control the onset and evolution of failure. These models are
used for delamination, adhesive failure, peeling, and other debonding
processes. They must be assigned within the material block of the
[`Interface` element](../elements/interface.md), which supplies the interface
kinematics and transfers the resulting normal and tangential tractions to the
adjacent continuum parts.

```{toctree}
:maxdepth: -1

dummy.md
powerlawmodei.md
thoulessmodei.md
xuneedleman.md
```
