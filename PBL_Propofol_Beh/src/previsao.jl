# =============================================================================
# previsao.jl — a previsão completa para um paciente, numa função só
#
# Junta as partes do modelo na ordem:
#   1. entrada escolhida pelo médico → taxa de infusão, segundo a segundo
#        propofol (mg/s): a) dose (mg/kg) em bolus + manutenção opcional (mg/kg/h); ou
#                         b) concentração-alvo no plasma (µg/ml)
#        remifentanil (µg/s), opcional: infusão contínua (µg/kg/min), começando
#                         um tempo antes do propofol, como nas induções do VitalDB
#   2. pacientes virtuais (populacao.jl) → uma curva de PAM por versão
#   3. faixa da PAM e períodos de risco (risco.jl)
# =============================================================================

"""
    taxa_infusao(idade, peso, altura, sexo; dose_mgkg, duracao_bolus_s, manutencao_mgkgh, alvo_ugml, duracao_s, opioide = false)

Taxa de infusão de propofol (mg/s), uma por segundo, de 0 a `duracao_s`, para a
entrada escolhida: dose (bolus + manutenção) ou concentração-alvo (a outra é
`nothing`). Na concentração-alvo, a bomba calcula com a farmacocinética típica
(com os termos de opioide se `opioide = true`).
"""
function taxa_infusao(idade, peso, altura, sexo; dose_mgkg, duracao_bolus_s, manutencao_mgkgh, alvo_ugml, duracao_s, opioide = false)
    amostras = duracao_s + 1
    if dose_mgkg !== nothing
        bolus = taxa_bolus(dose_mgkg * peso, duracao_bolus_s, amostras)
        # manutenção: mg/kg/h → mg/s (÷ 3600), só depois do bolus
        manut = [t < duracao_bolus_s ? 0.0 : manutencao_mgkgh * peso / 3600 for t in 0:duracao_s]
        return bolus .+ manut
    else
        return taxa_alvo(parametros_eleveld(idade, peso, altura, sexo; opioide), alvo_ugml, amostras)
    end
end

# Taxa do remifentanil (µg/s) de -antes_s até duracao_s, ou `nothing` se não houver remifentanil.
taxa_remi_ou_nada(peso, ugkgmin, antes_s, duracao_s) =
    ugkgmin > 0 ? taxa_remifentanil(peso; ugkgmin, antes_s, duracao_s) : nothing

"""
    prever(; idade, peso, altura, sexo, pam, fc,
             dose_mgkg = nothing, duracao_bolus_s = 30, manutencao_mgkgh = 0,
             alvo_ugml = nothing, remi_ugkgmin = 0, remi_antes_s = 540,
             duracao_s = 1200, n = 500, semente = 1)

Previsão da PAM durante a indução, como faixa, para um paciente.

Dados do paciente: `idade` (anos), `peso` (kg), `altura` (cm), `sexo` ("M" ou
"F"), `pam` e `fc` medidas antes de qualquer fármaco (mmHg e bpm).

Propofol (escolha UMA entrada):
- `dose_mgkg`: bolus em mg/kg, injetado em `duracao_bolus_s` segundos, seguido
  de `manutencao_mgkgh` mg/kg/h (0 = sem manutenção);
- `alvo_ugml`: concentração plasmática-alvo em µg/ml.

Remifentanil (opcional): `remi_ugkgmin` µg/kg/min em infusão contínua, começando
`remi_antes_s` segundos antes do propofol (0 = sem remifentanil). Padrão de 540 s:
nas 21 induções do VitalDB do PBL_6, o remifentanil começou uma mediana de 535 s
antes do propofol (dados/remifentanil_vitaldb.csv).

Devolve `(t, taxa, taxa_remi, faixa, periodos, texto, fc, cp, cr, rpt, vs, metricas, ...)`:
tempo (s, 0 = início do propofol), taxas de infusão, faixa da PAM (`inf`,
`mediana`, `sup`), períodos de risco e uma frase por período; faixas da FC, da
Cp do propofol e da Cp do remifentanil; medianas da RPT e do VS (fração do
basal); o resumo (`metricas`), o sono (BIS) e as farmacocinéticas típicas.
"""
function prever(; idade, peso, altura, sexo, pam, fc,
                  dose_mgkg = nothing, duracao_bolus_s = 30, manutencao_mgkgh = 0,
                  alvo_ugml = nothing, remi_ugkgmin = 0, remi_antes_s = 540,
                  duracao_s = 1200, n = 500, semente = 1)
    # conferências dos dados de entrada
    (dose_mgkg === nothing) == (alvo_ugml === nothing) &&
        error("Escolha UMA entrada: dose_mgkg (dose) ou alvo_ugml (concentração-alvo).")
    sexo in ("M", "F") || error("sexo deve ser \"M\" ou \"F\".")
    all(>(0), (idade, peso, altura, pam, fc)) || error("idade, peso, altura, pam e fc devem ser positivos.")
    remi_ugkgmin >= 0 && remi_antes_s >= 0 || error("O remifentanil e o tempo de início não podem ser negativos.")

    # 1. taxas de infusão: remifentanil (µg/s, de -remi_antes_s em diante) e propofol (mg/s, de 0 a duracao_s)
    taxa_r = taxa_remi_ou_nada(peso, remi_ugkgmin, remi_antes_s, duracao_s)
    taxa = taxa_infusao(idade, peso, altura, sexo; dose_mgkg, duracao_bolus_s, manutencao_mgkgh, alvo_ugml,
                        duracao_s, opioide = taxa_r !== nothing)

    # 2. pacientes virtuais e 3. faixa e períodos de risco
    pop = simular_populacao(idade, peso, altura, sexo, pam, fc, taxa; taxa_remi = taxa_r, antes_s = remi_antes_s, n, semente)
    faixa = faixa_pam(pop.curvas)
    periodos = periodos_de_risco(pop.t, faixa.inf, faixa.sup)
    # as outras saídas, também como faixa (FC e Cp) ou mediana (RPT e VS)
    faixa_fc, faixa_cp, faixa_cr = faixa_pam(pop.fc), faixa_pam(pop.cp), faixa_pam(pop.cr)
    rpt, vs = faixa_pam(pop.rpt).mediana, faixa_pam(pop.vs).mediana
    mt = metricas(pop.t, pam, fc, faixa, faixa_fc, faixa_cp, rpt, vs)
    # "metade dos pacientes": percentis 25% e 75% da PAM (a faixa de antes é de 5% a 95%)
    metade = faixa_pam(pop.curvas; inferior = 0.25, superior = 0.75)
    # sono: BIS como faixa (5–95% e 25–75%), % de pacientes dormindo (BIS < 60) a cada instante e,
    # pela mediana do BIS, quando o paciente "em média" dorme e acorda
    faixa_bis = faixa_pam(pop.bis)
    metade_bis = faixa_pam(pop.bis; inferior = 0.25, superior = 0.75)
    dormindo = [100 * count(b -> b[k] < LIMIAR_BIS, pop.bis) / length(pop.bis) for k in eachindex(pop.t)]
    sono = metricas_bis(pop.t, faixa_bis.mediana)
    i_acorda = sono.t_hipnose === nothing ? nothing : findnext(>=(LIMIAR_BIS), faixa_bis.mediana, findfirst(==(sono.t_hipnose), pop.t))
    sono = merge(sono, (t_acorda = i_acorda === nothing ? nothing : pop.t[i_acorda],
                        pct_max = maximum(dormindo), t_pct_max = pop.t[argmax(dormindo)]))
    return (t = pop.t, taxa = taxa, taxa_remi = taxa_r, faixa = faixa, periodos = periodos, texto = descrever(periodos),
            fc = faixa_fc, cp = faixa_cp, cr = faixa_cr, rpt = rpt, vs = vs, metricas = mt, metade = metade,
            bis = faixa_bis, metade_bis = metade_bis, dormindo = dormindo, sono = sono,
            pk = parametros_eleveld(idade, peso, altura, sexo; opioide = taxa_r !== nothing),
            pk_remi = parametros_minto(idade, peso, altura, sexo))
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
    simular_tipico(idade, peso, altura, sexo, pam, fc, taxa, taxa_remi, antes_s)

O paciente TÍPICO (os valores médios dos modelos, sem sorteio de versões): PAM
(mmHg) e BIS do início do propofol (0 s) até o fim. `taxa`: propofol (mg/s, de 0
em diante); `taxa_remi`: remifentanil (µg/s, de -antes_s em diante) ou `nothing`.
"""
function simular_tipico(idade, peso, altura, sexo, pam, fc, taxa, taxa_remi, antes_s)
    com_remi = taxa_remi !== nothing
    antes = com_remi ? antes_s : 0
    pk = parametros_eleveld(idade, peso, altura, sexo; opioide = com_remi)
    tp = vcat(zeros(antes), taxa)
    cr = com_remi ? simular_pk(parametros_minto(idade, peso, altura, sexo), taxa_remi) : nothing
    pam_t = simular_hemodinamica(simular_pk(pk, tp), pam, fc, idade; cr, subpassos = 1).pam
    bis = simular_bis(pk, tp, idade, peso).bis
    return (pam = pam_t[antes + 1:end], bis = bis[antes + 1:end])
end

# Uma linha das tabelas de pressão × dose: PAM mínima, quando, queda, segundos abaixo
# de 65 e de 55, nível de risco (PAM mínima arredondada, como no monitor) e BIS mínimo.
function linha_tipico(valor, total, r, pam; limiares = LIMIARES)
    i = argmin(r.pam)
    return (valor = valor, total = total, pam_min = r.pam[i], t_min = i - 1.0, queda_pct = 100 * (1 - r.pam[i] / pam),
            s_abaixo_65 = count(<(limiares.moderado), r.pam) * 1.0, s_abaixo_55 = count(<(limiares.alto), r.pam) * 1.0,
            risco = nivel(round(r.pam[i]); limiares), bis_min = minimum(r.bis))
end

"""
    risco_por_dose(valores; idade, peso, altura, sexo, pam, fc, modo = "dose",
                   duracao_bolus_s = 30, manutencao_mgkgh = 0, remi_ugkgmin = 0, remi_antes_s = 540,
                   duracao_s = 1200, limiares = LIMIARES)

Gráfico pressão × dose de PROPOFOL, com o remifentanil fixo (`remi_ugkgmin`): a
PAM mínima do paciente típico para cada valor de `valores` — doses em mg/kg
(`modo = "dose"`) ou concentrações-alvo em µg/ml (`modo = "alvo"`). Uma simulação
só por valor: é rápida e permite varrer muitas doses de uma vez.

Devolve, para cada valor: `valor`, `total` (mg de propofol), `pam_min` (mmHg),
`t_min` (s), `queda_pct`, `s_abaixo_65` e `s_abaixo_55` (s), `risco` e `bis_min`.
"""
function risco_por_dose(valores; idade, peso, altura, sexo, pam, fc, modo = "dose",
                        duracao_bolus_s = 30, manutencao_mgkgh = 0, remi_ugkgmin = 0, remi_antes_s = 540,
                        duracao_s = 1200, limiares = LIMIARES)
    taxa_r = taxa_remi_ou_nada(peso, remi_ugkgmin, remi_antes_s, duracao_s)
    map(valores) do v
        taxa = taxa_infusao(idade, peso, altura, sexo; dose_mgkg = modo == "alvo" ? nothing : v, duracao_bolus_s,
                            manutencao_mgkgh, alvo_ugml = modo == "alvo" ? v : nothing, duracao_s, opioide = taxa_r !== nothing)
        linha_tipico(v, sum(taxa), simular_tipico(idade, peso, altura, sexo, pam, fc, taxa, taxa_r, remi_antes_s), pam; limiares)
    end
end

"""
    risco_por_remi(valores; idade, peso, altura, sexo, pam, fc, modo = "dose", dose_mgkg = nothing,
                   alvo_ugml = nothing, duracao_bolus_s = 30, manutencao_mgkgh = 0, remi_antes_s = 540,
                   duracao_s = 1200, limiares = LIMIARES)

Gráfico pressão × dose de REMIFENTANIL, com o propofol fixo (`dose_mgkg` ou
`alvo_ugml`): a PAM mínima do paciente típico para cada taxa de remifentanil em
`valores` (µg/kg/min). Devolve as mesmas colunas de `risco_por_dose`, com `total`
em µg de remifentanil.
"""
function risco_por_remi(valores; idade, peso, altura, sexo, pam, fc, dose_mgkg = nothing, alvo_ugml = nothing,
                        duracao_bolus_s = 30, manutencao_mgkgh = 0, remi_antes_s = 540, duracao_s = 1200, limiares = LIMIARES)
    map(valores) do v
        taxa_r = taxa_remi_ou_nada(peso, v, remi_antes_s, duracao_s)
        taxa = taxa_infusao(idade, peso, altura, sexo; dose_mgkg, duracao_bolus_s, manutencao_mgkgh, alvo_ugml,
                            duracao_s, opioide = taxa_r !== nothing)
        total = taxa_r === nothing ? 0.0 : sum(taxa_r)
        linha_tipico(v, total, simular_tipico(idade, peso, altura, sexo, pam, fc, taxa, taxa_r, remi_antes_s), pam; limiares)
    end
end

"""
    chances_por_dose(valores; idade, peso, altura, sexo, pam, fc, modo = "dose",
                     duracao_bolus_s = 30, manutencao_mgkgh = 0, remi_ugkgmin = 0, remi_antes_s = 540,
                     duracao_s = 1200, n = 300, semente = 1)

Para o gráfico de resposta à dose: para cada valor (dose em mg/kg ou alvo em µg/ml),
com o remifentanil fixo, os mesmos `n` pacientes virtuais (mesma `semente`) e a
porcentagem deles que dorme (BIS < 60 em algum momento), que tem PAM < 65 e < 55
mmHg (PAM arredondada, como no monitor) e que dorme sem PAM < 65 (`ok`).
"""
function chances_por_dose(valores; idade, peso, altura, sexo, pam, fc, modo = "dose",
                          duracao_bolus_s = 30, manutencao_mgkgh = 0, remi_ugkgmin = 0, remi_antes_s = 540,
                          duracao_s = 1200, n = 300, semente = 1)
    taxa_r = taxa_remi_ou_nada(peso, remi_ugkgmin, remi_antes_s, duracao_s)
    map(valores) do v
        taxa = taxa_infusao(idade, peso, altura, sexo; dose_mgkg = modo == "alvo" ? nothing : v, duracao_bolus_s,
                            manutencao_mgkgh, alvo_ugml = modo == "alvo" ? v : nothing, duracao_s, opioide = taxa_r !== nothing)
        pop = simular_populacao(idade, peso, altura, sexo, pam, fc, taxa; taxa_remi = taxa_r, antes_s = remi_antes_s, n, semente)
        dorme = [minimum(b) < LIMIAR_BIS for b in pop.bis]
        pmin = [round(minimum(c)) for c in pop.curvas]
        pct(x) = 100 * count(x) / n
        (valor = v, dorme = pct(dorme), h65 = pct(pmin .< 65), h55 = pct(pmin .< 55), ok = pct(dorme .& (pmin .>= 65)))
    end
end
