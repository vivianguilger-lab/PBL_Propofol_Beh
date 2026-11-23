# =============================================================================
# remifentanil.jl — da infusão de remifentanil à concentração no sangue (Cp)
#
# Modelo: Minto CF, Schnider TW, Egan TD et al. Influence of age and gender on
# the pharmacokinetics and pharmacodynamics of remifentanil: I. Model development.
# Anesthesiology 1997;86(1):10–23. doi:10.1097/00000542-199701000-00004
#
# Por que Minto:
#   - é o modelo que a Beatriz implementou no PBL_6 (src/remifentanil.jl, mesmas
#     equações) e o padrão do PAS para o remifentanil (pk_models.py L357–L371);
#   - reproduz a Cp que as bombas do VitalDB mostraram: com a taxa registrada
#     (Orchestra/RFTN20_RATE), a Cp de Minto difere da Cp da bomba
#     (Orchestra/RFTN20_CP) em 1% (mediana) nos 21 casos da coorte do PBL_6
#     (scripts/remifentanil_vitaldb.py, dados/remifentanil_vitaldb.csv).
#
# Unidades: tempo em segundos; dose em µg; volumes em L; Cp em µg/L = ng/ml.
# A conta é a mesma do propofol (simular_pk, farmacocinetica.jl): 3 compartimentos.
# =============================================================================

# Massa magra de James, a que o modelo de Minto usa (PAS patient.py L217–L219).
massa_magra_james(peso, altura, sexo) =
    sexo == "M" ? 1.1 * peso - 128 * (peso / altura)^2 : 1.07 * peso - 148 * (peso / altura)^2

"""
    parametros_minto(idade, peso, altura, sexo)

Volumes (L) e depurações (L/min) do modelo de Minto para um paciente
(PAS pk_models.py L362–L371; PBL_6 src/remifentanil.jl). Referência do modelo:
40 anos e 55 kg de massa magra.
"""
function parametros_minto(idade, peso, altura, sexo)
    Δidade = idade - 40
    Δmagra = massa_magra_james(peso, altura, sexo) - 55
    pk = (V1 = 5.1 - 0.0201 * Δidade + 0.072 * Δmagra,     # sangue (L)
          V2 = 9.82 - 0.0811 * Δidade + 0.108 * Δmagra,    # compartimento rápido (L)
          V3 = 5.42,                                       # compartimento lento (L)
          CL = 2.6 - 0.0162 * Δidade + 0.0191 * Δmagra,    # eliminação (L/min)
          Q2 = 2.05 - 0.0301 * Δidade,                     # troca sangue ⇄ rápido (L/min)
          Q3 = 0.076 - 0.00113 * Δidade)                   # troca sangue ⇄ lento (L/min)
    all(>(0), pk) || error("Os dados do paciente estão fora da faixa do modelo de Minto para o remifentanil.")
    return pk
end

# Variabilidade entre pessoas (ω², log-normal) do modelo de Minto: o PAS guarda os
# coeficientes de variação (pk_models.py L375–L380) e converte em ω² = ln(1 + CV²) (L384–L389).
const ω²_MINTO = (V1 = log(1 + 0.26^2), V2 = log(1 + 0.29^2), V3 = log(1 + 0.66^2),
                  CL = log(1 + 0.14^2), Q2 = log(1 + 0.36^2), Q3 = log(1 + 0.41^2))

"""
    taxa_remifentanil(peso; ugkgmin, antes_s, duracao_s)

Infusão contínua de remifentanil (µg/s), uma por segundo, de `-antes_s` (início
do remifentanil, antes do propofol) até `duracao_s` (fim da simulação). O propofol
começa em 0. `ugkgmin`: taxa em µg/kg/min, como na bomba.
"""
taxa_remifentanil(peso; ugkgmin, antes_s, duracao_s) = fill(ugkgmin * peso / 60, antes_s + duracao_s + 1)
