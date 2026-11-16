# Gera valores de referência com o PAS (commit 4afc03a): Cp do propofol (Eleveld,
# do próprio PAS) e a resposta hemodinâmica de Su 2023 com remifentanil = 0.
# Uso (precisa de numpy, scipy, matplotlib e casadi):
#   git clone https://github.com/AnesthesiaSimulation/Python_Anesthesia_Simulator pas
#   git -C pas checkout 4afc03ad2997f725205c202e0d0634ff0664e5c7
#   python gerar_referencia_pas.py pas/src referencia_pas.csv
import sys, csv
import numpy as np
from scipy.integrate import solve_ivp
sys.path.insert(0, sys.argv[1])
from python_anesthesia_simulator.pk_models import CompartmentModel
from python_anesthesia_simulator.pd_models import Hemo_meca_PD_model

idade, altura, peso, sexo = 62, 160, 72, 0          # sexo 0 = feminino no PAS
lbm = 1.07 * peso - 148 * (peso / altura) ** 2       # James (mulher), só para o PK do PAS
pk = CompartmentModel([idade, altura, peso, sexo], lbm, drug="Propofol", model="Eleveld", ts=1)
dose_mg = 1.5 * peso
taxa = dose_mg / 30                                   # mg/s durante 30 s
cp = []
for k in range(1200):
    cp.append(float(pk.x[0, 0]))
    pk.one_step(taxa if k < 30 else 0.0)
cp = np.array(cp)

def roda_pas(sv):
    h = Hemo_meca_PD_model(age=idade, ts=1, hr_base=78, sv_base=sv, map_base=104)
    out = [h.output_function(h.x)]
    for k in range(1, len(cp)):
        out.append(h.one_step(cp_propo=cp[k], cp_remi=0.0))
    return np.array(out)

def roda_exato(sv):
    # mesmas equações do PAS, integradas com tolerância bem apertada
    h = Hemo_meca_PD_model(age=idade, ts=1, hr_base=78, sv_base=sv, map_base=104)
    x = h.x.copy(); out = [h.output_function(x)]
    for k in range(1, len(cp)):
        u = [(cp[k - 1] + cp[k]) / 2, 0.0, 0, 0, 0, 0, 0]
        r = solve_ivp(h.continuous_dynamic_sys, [0, 1], x, args=(u,), rtol=1e-11, atol=1e-13)
        x = r.y[:, -1]; out.append(h.output_function(x))
    return np.array(out)

a, b, e = roda_pas(70.0), roda_pas(90.0), roda_exato(70.0)
print("PAM mínima PAS (VS 70): %.4f | (VS 90): %.4f | exato: %.4f" % (a[:, 3].min(), b[:, 3].min(), e[:, 3].min()))
print("maior diferença PAS×exato: PAM %.2e  FC %.2e" % (abs(a[:, 3] - e[:, 3]).max(), abs(a[:, 2] - e[:, 2]).max()))
with open(sys.argv[2], "w", newline="") as f:
    w = csv.writer(f); w.writerow(["t_s", "cp", "pam_pas", "fc_pas", "pam_exato", "fc_exato", "vs_rel_exato", "rpt_rel_exato"])
    for k in range(len(cp)):
        w.writerow([k, cp[k], a[k, 3], a[k, 2], e[k, 3], e[k, 2], e[k, 1] / e[0, 1], e[k, 0] / e[0, 0]])
