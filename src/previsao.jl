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
    taxa_infusao(idade, peso, altura, sexo; dose_mgkg, duracao_bolus_s, manutencao_mgkgh, alvo_ugml, duracao_s)

Taxa de infusão (mg/s), uma por segundo, de 0 a `duracao_s`, para a entrada
escolhida: dose (bolus + manutenção) ou concentração-alvo (a outra é `nothing`).
"""
function taxa_infusao(idade, peso, altura, sexo; dose_mgkg, duracao_bolus_s, manutencao_mgkgh, alvo_ugml, duracao_s)
    amostras = duracao_s + 1
    if dose_mgkg !== nothing
        bolus = taxa_bolus(dose_mgkg * peso, duracao_bolus_s, amostras)
        # manutenção: mg/kg/h → mg/s (÷ 3600), só depois do bolus
        manut = [t < duracao_bolus_s ? 0.0 : manutencao_mgkgh * peso / 3600 for t in 0:duracao_s]
        return bolus .+ manut
    else
        return taxa_alvo(parametros_eleveld(idade, peso, altura, sexo), alvo_ugml, amostras)
    end
end

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

Devolve `(t, taxa, faixa, periodos, texto, fc, cp, rpt, vs, metricas, pk)`:
tempo (s), taxa de infusão (mg/s), faixa da PAM (`inf`, `mediana`, `sup`),
períodos de risco e uma frase por período; faixas da FC e da Cp; medianas da RPT
e do VS (fração do basal); o resumo (`metricas`) e a farmacocinética típica (`pk`).
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
    taxa = taxa_infusao(idade, peso, altura, sexo; dose_mgkg, duracao_bolus_s, manutencao_mgkgh, alvo_ugml, duracao_s)

    # 2. pacientes virtuais e 3. faixa e períodos de risco
    pop = simular_populacao(idade, peso, altura, sexo, pam, fc, taxa; n, semente)
    faixa = faixa_pam(pop.curvas)
    periodos = periodos_de_risco(pop.t, faixa.inf, faixa.sup)
    # as outras saídas, também como faixa (FC e Cp) ou mediana (RPT e VS)
    faixa_fc, faixa_cp = faixa_pam(pop.fc), faixa_pam(pop.cp)
    rpt, vs = faixa_pam(pop.rpt).mediana, faixa_pam(pop.vs).mediana
    mt = metricas(pop.t, pam, fc, faixa, faixa_fc, faixa_cp, rpt, vs)
    return (t = pop.t, taxa = taxa, faixa = faixa, periodos = periodos, texto = descrever(periodos),
            fc = faixa_fc, cp = faixa_cp, rpt = rpt, vs = vs, metricas = mt,
            pk = parametros_eleveld(idade, peso, altura, sexo))
end

"""
    metricas(t, pam_basal, fc_basal, faixa, fc, cp, rpt, vs; limiares = LIMIARES)

Resumo da previsão, para o cartão de risco e os blocos de números:
- PAM mínima da mediana (`pam_nadir`, em `t_nadir` s) e queda em % da basal;
- pior caso: PAM mínima da borda de baixo da faixa (`pam_pior`, em `t_pior` s) e o
  nível de risco dele (`risco`, com a PAM arredondada, como em `periodos_de_risco`);
- segundos com PAM abaixo de 65 e de 55, na mediana e no pior caso;
- FC mínima e Cp máxima (medianas);
- variação (%) de RPT, VS e FC no instante da PAM mínima da mediana.
"""
function metricas(t, pam_basal, fc_basal, faixa, fc, cp, rpt, vs; limiares = LIMIARES)
    dt = t[2] - t[1]
    abaixo(curva, l) = dt * count(<(l), curva)        # segundos abaixo do limiar
    i, j, k = argmin(faixa.mediana), argmin(faixa.inf), argmax(cp.mediana)
    return (pam_basal = pam_basal, pam_nadir = faixa.mediana[i], t_nadir = t[i],
            queda_pct = 100 * (1 - faixa.mediana[i] / pam_basal),
            pam_pior = faixa.inf[j], t_pior = t[j], risco = nivel(round(faixa.inf[j]); limiares),
            s_abaixo_65 = abaixo(faixa.mediana, limiares.moderado), s_abaixo_65_pior = abaixo(faixa.inf, limiares.moderado),
            s_abaixo_55 = abaixo(faixa.mediana, limiares.alto), s_abaixo_55_pior = abaixo(faixa.inf, limiares.alto),
            fc_min = minimum(fc.mediana), cp_max = cp.mediana[k], t_cp_max = t[k],
            rpt_nadir_pct = 100 * (rpt[i] - 1), vs_nadir_pct = 100 * (vs[i] - 1),
            fc_nadir_pct = 100 * (fc.mediana[i] / fc_basal - 1))
end

"""
    risco_por_dose(valores; idade, peso, altura, sexo, pam, fc, modo = "dose",
                   duracao_bolus_s = 30, manutencao_mgkgh = 0, duracao_s = 1200, limiares = LIMIARES)

Para o gráfico pressão × dose: a PAM mínima do paciente TÍPICO (sem sorteio de
versões) para cada valor de `valores` — doses em mg/kg (`modo = "dose"`) ou
concentrações-alvo em µg/ml (`modo = "alvo"`). Como é uma simulação só por
valor, é rápida e permite varrer muitas doses de uma vez.

Devolve, para cada valor: `valor`, `mg` (total infundido), `pam_min` (mmHg),
`t_min` (s), `queda_pct`, `s_abaixo_65` e `s_abaixo_55` (s) e `risco` (nível da
PAM mínima arredondada).
"""
function risco_por_dose(valores; idade, peso, altura, sexo, pam, fc, modo = "dose",
                        duracao_bolus_s = 30, manutencao_mgkgh = 0, duracao_s = 1200, limiares = LIMIARES)
    pk = parametros_eleveld(idade, peso, altura, sexo)
    map(valores) do v
        taxa = modo == "alvo" ?
            taxa_infusao(idade, peso, altura, sexo; dose_mgkg = nothing, duracao_bolus_s, manutencao_mgkgh, alvo_ugml = v, duracao_s) :
            taxa_infusao(idade, peso, altura, sexo; dose_mgkg = v, duracao_bolus_s, manutencao_mgkgh, alvo_ugml = nothing, duracao_s)
        curva = simular_hemodinamica(simular_pk(pk, taxa), pam, fc, idade).pam
        i = argmin(curva)
        (valor = v, mg = sum(taxa), pam_min = curva[i], t_min = i - 1.0, queda_pct = 100 * (1 - curva[i] / pam),
         s_abaixo_65 = count(<(limiares.moderado), curva) * 1.0, s_abaixo_55 = count(<(limiares.alto), curva) * 1.0,
         risco = nivel(round(curva[i]); limiares))
    end
end
