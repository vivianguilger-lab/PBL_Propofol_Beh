# Gera valores de referência com o PAS (commit 4afc03a) para o remifentanil:
# Cp do remifentanil (Minto, do próprio PAS), Cp do propofol (Eleveld com opioide)
# e a resposta hemodinâmica de Su 2023 com propofol E remifentanil.
# Cenário: mulher, 62 anos, 72 kg, 160 cm, PAM 104 e FC 78 antes de qualquer fármaco;
# remifentanil 0,12 µg/kg/min começando 540 s antes do propofol; propofol 1,5 mg/kg em 30 s.
# Uso (precisa de numpy, scipy, matplotlib e casadi):
#   git clone https://github.com/AnesthesiaSimulation/Python_Anesthesia_Simulator pas
#   git -C pas checkout 4afc03ad2997f725205c202e0d0634ff0664e5c7
#   python gerar_referencia_remi.py pas/src referencia_remi.csv
import sys, csv
import numpy as np
from scipy.integrate import solve_ivp
sys.path.insert(0, sys.argv[1])
from python_anesthesia_simulator.pk_models import CompartmentModel
from python_anesthesia_simulator.pd_models import Hemo_meca_PD_model

idade, altura, peso, sexo = 62, 160, 72, 0          # sexo 0 = feminino no PAS
lbm = 1.07 * peso - 148 * (peso / altura) ** 2       # James (mulher): o PAS usa para o Minto
antes, dur = 540, 1200
pk_p = CompartmentModel([idade, altura, peso, sexo], lbm, drug="Propofol", model="Eleveld", ts=1, opiate=True)
pk_r = CompartmentModel([idade, altura, peso, sexo], lbm, drug="Remifentanil", model="Minto", ts=1)
taxa_p = 1.5 * peso / 30                              # mg/s durante 30 s
taxa_r = 0.12 * peso / 60                             # µg/s, contínua
cp, cr = [], []
for k in range(antes + dur + 1):
    cp.append(float(pk_p.x[0, 0])); cr.append(float(pk_r.x[0, 0]))
    t = k - antes                                     # t = 0: início do propofol
    pk_p.one_step(taxa_p if 0 <= t < 30 else 0.0)
    pk_r.one_step(taxa_r)
cp, cr = np.array(cp), np.array(cr)

# mesmas equações do PAS, integradas com tolerância bem apertada
h = Hemo_meca_PD_model(age=idade, ts=1, hr_base=78, sv_base=70, map_base=104)
x = h.x.copy(); out = [h.output_function(x)]
for k in range(1, len(cp)):
    u = [(cp[k - 1] + cp[k]) / 2, (cr[k - 1] + cr[k]) / 2, 0, 0, 0, 0, 0]
    r = solve_ivp(h.continuous_dynamic_sys, [0, 1], x, args=(u,), rtol=1e-11, atol=1e-13)
    x = r.y[:, -1]; out.append(h.output_function(x))
e = np.array(out)
print("PAM no início do propofol: %.3f | mínima: %.3f | Cp remi no início do propofol: %.3f ng/ml"
      % (e[antes, 3], e[:, 3].min(), cr[antes]))
with open(sys.argv[2], "w", newline="") as f:
    w = csv.writer(f); w.writerow(["t_s", "cp", "cr", "pam", "fc"])
    for k in range(len(cp)):
        w.writerow([k - antes, cp[k], cr[k], e[k, 3], e[k, 2]])
