# Gera valores de referência do BIS de Eleveld com o PAS (commit 4afc03a): concentração
# no sítio efetor (Ce, µg/ml) e BIS com o atraso do monitor, segundo a segundo, sem opioide.
# Uso (precisa de numpy, scipy, matplotlib e casadi):
#   git clone https://github.com/AnesthesiaSimulation/Python_Anesthesia_Simulator pas
#   git -C pas checkout 4afc03ad2997f725205c202e0d0634ff0664e5c7
#   python gerar_referencia_bis.py pas/src referencia_bis.csv
import sys, csv
sys.path.insert(0, sys.argv[1])
from python_anesthesia_simulator.pk_models import CompartmentModel
from python_anesthesia_simulator.pd_models import BIS_model

N = 1200  # segundos
# Paciente A: mulher, 62 anos, 72 kg, 160 cm; 1,5 mg/kg em 30 s.
# Paciente B: homem, 30 anos, 80 kg, 180 cm; 2 mg/kg em 10 s e depois 100 µg/kg/min.
pacientes = {
    "A": dict(idade=62, altura=160, peso=72, sexo=0,
              taxa=[1.5 * 72 / 30 if t < 30 else 0.0 for t in range(N)]),
    "B": dict(idade=30, altura=180, peso=80, sexo=1,
              taxa=[2 * 80 / 10 if t < 10 else 0.1 * 80 / 60 for t in range(N)]),
}
colunas = {}
for nome, p in pacientes.items():
    pk = CompartmentModel([p["idade"], p["altura"], p["peso"], p["sexo"]], 50.0,
                          drug="Propofol", model="Eleveld", ts=1, opiate=False)
    bis_pd = BIS_model(hill_model="Eleveld", age=p["idade"], ts=1)
    ce, bis = [], []
    for k in range(N):
        c = float(pk.x[3, 0])          # Ce no início do segundo k (4º compartimento: sítio efetor)
        ce.append(c)
        bis.append(float(bis_pd.one_step(c)))   # o PAS devolve o BIS com o atraso do monitor
        pk.one_step(p["taxa"][k])
    colunas["taxa_" + nome] = p["taxa"]
    colunas["ce_" + nome] = ce
    colunas["bis_" + nome] = bis
with open(sys.argv[2], "w", newline="") as f:
    w = csv.writer(f)
    nomes = list(colunas)
    w.writerow(["t_s"] + nomes)
    for k in range(N):
        w.writerow([k] + [colunas[c][k] for c in nomes])
