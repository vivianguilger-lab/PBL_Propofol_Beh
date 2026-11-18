# =============================================================================
# previsao.jl — a previsão completa para um paciente, numa função só
#
# Junta as partes do modelo na ordem:
#   1. entrada escolhida pelo médico → taxa de infusão (mg/s), segundo a segundo
#        a) dose (mg/kg) em bolus + manutenção opcional (mg/kg/h); ou
#        b) concentração-alvo no plasma (µg/ml)
#   2. pacientes virtuais (populacao.jl) → uma curva de PAM por versão
#   3. faixa da PAM e períodos de risco (risco.jl)
# =============================================================================

"""
    prever(; idade, peso, altura, sexo, pam, fc,
             dose_mgkg = nothing, duracao_bolus_s = 30, manutencao_mgkgh = 0,
             alvo_ugml = nothing, duracao_s = 1200, n = 500, semente = 1)

Previsão da PAM durante a indução, como faixa, para um paciente.

Dados do paciente: `idade` (anos), `peso` (kg), `altura` (cm), `sexo` ("M" ou
"F"), `pam` e `fc` medidas antes da indução (mmHg e bpm).

Entrada (escolha UMA):
- `dose_mgkg`: bolus em mg/kg, injetado em `duracao_bolus_s` segundos, seguido
  de `manutencao_mgkgh` mg/kg/h (0 = sem manutenção);
- `alvo_ugml`: concentração plasmática-alvo em µg/ml.

Devolve `(t, taxa, faixa, periodos, texto)`: tempo (s), taxa de infusão (mg/s),
faixa da PAM (`inf`, `mediana`, `sup`), períodos de risco e uma frase por período.
"""
function prever(; idade, peso, altura, sexo, pam, fc,
                  dose_mgkg = nothing, duracao_bolus_s = 30, manutencao_mgkgh = 0,
                  alvo_ugml = nothing, duracao_s = 1200, n = 500, semente = 1)
    # conferências dos dados de entrada
    (dose_mgkg === nothing) == (alvo_ugml === nothing) &&
        error("Escolha UMA entrada: dose_mgkg (dose) ou alvo_ugml (concentração-alvo).")
    sexo in ("M", "F") || error("sexo deve ser \"M\" ou \"F\".")
    all(>(0), (idade, peso, altura, pam, fc)) || error("idade, peso, altura, pam e fc devem ser positivos.")

    # 1. taxa de infusão (mg/s), uma por segundo, de 0 a duracao_s
    amostras = duracao_s + 1
    if dose_mgkg !== nothing
        bolus = taxa_bolus(dose_mgkg * peso, duracao_bolus_s, amostras)
        # manutenção: mg/kg/h → mg/s (÷ 3600), só depois do bolus
        manut = [t < duracao_bolus_s ? 0.0 : manutencao_mgkgh * peso / 3600 for t in 0:duracao_s]
        taxa = bolus .+ manut
    else
        taxa = taxa_alvo(parametros_eleveld(idade, peso, altura, sexo), alvo_ugml, amostras)
    end

    # 2. pacientes virtuais e 3. faixa e períodos de risco
    pop = simular_populacao(idade, peso, altura, sexo, pam, fc, taxa; n, semente)
    faixa = faixa_pam(pop.curvas)
    periodos = periodos_de_risco(pop.t, faixa.inf, faixa.sup)
    return (t = pop.t, taxa = taxa, faixa = faixa, periodos = periodos, texto = descrever(periodos))
end
