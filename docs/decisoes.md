# Decisões do projeto

Cada decisão de método, com o motivo, a fonte e a situação. **Pendente** quer
dizer que ainda falta conferir ou decidir.

## Escopo

| # | Decisão | Por quê / fonte | Situação |
|---|---|---|---|
| 1 | Só propofol. O remifentanil entra com **concentração zero** (os coeficientes do Su ficam como publicados). | Foco do projeto. Com concentração 0, a curva do remifentanil, C/(C50 + C), dá 0: não há divisão por zero e o efeito do propofol fica igual ao do artigo. Zerar os *coeficientes* quebraria o modelo (C50 = 0 dá 0/0). | **Pendente:** dúvida sobre usar dados de indução combinada do propofol e remifetanil ou tentar achar só do propofol. Conferir no artigo se houve períodos só com propofol. |

## Modelo

| # | Decisão | Por quê / fonte | Situação |
|---|---|---|---|
| 2 | Efeito do propofol na circulação: estrutura e valores de Su 2023. | Su H et al. *Br J Anaesth* 2023;131:222–33; implementação de referência no PAS (`pd_models.py` L890). Valores estimados em voluntários saudáveis: ponto de partida, não verdade final. | Verificado contra o PAS (diferença < 10⁻¹² mmHg). **Pendente:** conferir os 12 valores nas tabelas do artigo; o PAS declara `age_max_sv = 0.033` mas usa 0,0333; significado do "termo dependente do tempo". |
| 3 | O VS basal não é usado. | A PAM e a FC previstas não dependem dele (tudo em fração do basal); no PAS, VS basal de 70 ou 90 ml dá a mesma PAM. | Confirmado (teste). |
| 4 | Farmacocinética: Eleveld 2018, concentração arterial, sem opioide. | Eleveld DJ et al. *Br J Anaesth* 2018;120:942–59. Mesmo grupo do Su; usado junto com o Su no PAS; população ampla; publica a variabilidade. | Verificado contra o PAS (diferença < 10⁻¹³ µg/ml). **Pendente:** qual PK o Su usou para calcular a Cp. |
| 5 | Entrada por concentração-alvo: alvo no **plasma**; taxa calculada com o paciente **típico**; sem limite de bomba por padrão. | O Su usa a Cp, e o alvo no sítio efetor exigiria um ke0. Uma bomba real calcula com o paciente típico. O limite da bomba depende do equipamento. | Testado (`taxa_alvo`). |
| 6 | Tempo em segundos em todo o código e nos gráficos. | Pedido do grupo, para facilitar a análise. | — |

## Faixas e risco

| # | Decisão | Por quê / fonte | Situação |
|---|---|---|---|
| 7 | Variabilidade entre pessoas: ω² de V1, V2, V3, CL, Q2 e Q3 (Eleveld) e do C50 da RPT (Su). Os outros parâmetros do Su ficam fixos. | É a variabilidade publicada que o PAS traz (`pk_models.py` L346–L351; `pd_models.py` L1009). A faixa provavelmente é **mais estreita** que a realidade. | **Pendente:** comparar com pacientes reais. |
| 8 | Sorteio do C50 da RPT pela definição do artigo (η com média 0 e variância ω²), e não pelo código do PAS. | O PAS passa ω como **média** do sorteio (`pd_models.py` L1054 e L1075); na farmacocinética (L518) ele faz certo. | Testado (média e variância). |
| 9 | Faixa = percentis 5% a 95% (90% das versões). | Convenção, não valor de artigo. Ajustável (`inferior`, `superior`). | — |
| 10 | Limiares de risco: baixo ≥ 65 mmHg; moderado de 55 a 65; alto < 55. Bordas da faixa arredondadas para mmHg inteiros antes de classificar. | 65: Salmasi V et al. *Anesthesiology* 2017;126:47–65. 55: Walsh M et al. *Anesthesiology* 2013;119:507–15. O arredondamento segue o monitor, que mostra mmHg inteiros, e evita períodos picotados de 1 a 4 s. | — |

## Dados reais (VitalDB)

| # | Decisão | Por quê / fonte | Situação |
|---|---|---|---|
| 11 | A seleção e o processamento dos casos do VitalDB ficam com a Beatriz, no PBL_6. Este projeto só vai **ler** os arquivos gerados lá. | Divisão de trabalho do grupo. | **Pendente (dúvida):** os 21 casos selecionados estão todos com `coadmin_status = active` (remifentanil junto). Usamos esses casos incluindo o remifentanil no modelo, ou buscamos casos sem remifentanil na indução? |

Contagens feitas no VitalDB (6388 casos), para ajudar nessa decisão:
- 3512 casos com bomba de propofol; **3509** deles também com remifentanil. Numa amostra de 76 casos, o remifentanil começou sempre antes ou junto com o propofol (mediana: 30 s antes).
- 80 casos de anestesia inalatória com propofol em bolus manual, sem fentanil nem midazolam. O horário do bolus não é registrado (só a dose total) e só 6 têm pressão invasiva antes da indução.
