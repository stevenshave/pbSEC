from high_accuracy_binding_equations import one_to_one_binding
cpd=8
kd=10
p=10
loss=1.0

while cpd>0.012:
    pl=one_to_one_binding(p,cpd,kd)
    cpd=float(pl.real*loss)
    p=p*loss
    print(f"{cpd=}, {p=}")