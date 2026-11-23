# Decisões do projeto

Cada decisão de método, com o motivo, a fonte e a situação. **Pendente** quer
dizer que ainda falta conferir ou decidir.

## Escopo

| # | Decisão | Por quê / fonte | Situação |
|---|---|---|---|
| 1 | ~~Só propofol. O remifentanil entra com concentração zero.~~ **Substituída pelos itens 12 a 19:** o remifentanil agora é escolhido pelo médico (0 = só propofol). | Nas 21 induções do VitalDB selecionadas no PBL_6, todas têm remifentanil junto, começando antes do propofol (item 18). Com remifentanil 0, o resultado é exatamente o de antes (teste). | Resolvida. |

## Modelo

| # | Decisão | Por quê / fonte | Situação |
|---|---|---|---|
| 2 | Efeito do propofol na circulação: estrutura e valores de Su 2023. | Su H et al. *Br J Anaesth* 2023;131:222–33; implementação de referência no PAS (`pd_models.py` L890). Valores estimados em voluntários saudáveis: ponto de partida, não verdade final. | Verificado contra o PAS (diferença < 10⁻¹² mmHg). **Pendente:** conferir os 12 valores nas tabelas do artigo; o PAS declara `age_max_sv = 0.033` mas usa 0,0333; significado do "termo dependente do tempo". |
| 3 | O VS basal não é usado. | A PAM e a FC previstas não dependem dele (tudo em fração do basal); no PAS, VS basal de 70 ou 90 ml dá a mesma PAM. | Confirmado (teste). |
| 4 | Farmacocinética: Eleveld 2018, concentração arterial; **com os termos de opioide quando há remifentanil** (item 15). | Eleveld DJ et al. *Br J Anaesth* 2018;120:942–59. Mesmo grupo do Su; usado junto com o Su no PAS; população ampla; publica a variabilidade. | Verificado contra o PAS (diferença < 10⁻¹³ µg/ml). **Pendente:** qual PK o Su usou para calcular a Cp. |
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

## Remifentanil

| # | Decisão | Por quê / fonte | Situação |
|---|---|---|---|
| 12 | Farmacocinética do remifentanil: **Minto 1997** (massa magra de James). | Minto CF et al. *Anesthesiology* 1997;86:10–23. É o modelo do `src/remifentanil.jl` da Beatriz no PBL_6 e o padrão do PAS (`pk_models.py` L357–L371). Com as taxas registradas no VitalDB (`Orchestra/RFTN20_RATE`), a Cp de Minto difere da Cp mostrada pela própria bomba (`Orchestra/RFTN20_CP`) em **0,97%** (mediana; de 0,21% a 4,0%) nos 21 casos (`scripts/remifentanil_vitaldb.py`, `dados/remifentanil_vitaldb.csv`). | Verificado contra o PAS (diferença < 10⁻¹³ ng/ml) e contra a bomba do VitalDB. **Pendente:** qual PK do remifentanil o Su usou (Minto ou Eleveld 2017). |
| 13 | Efeito do remifentanil na circulação: os 8 coeficientes do remifentanil e das interações de **Su 2023**, como publicados. | Su 2023; PAS `pd_models.py` L971–L978 (valores) e L1149–L1182 (equações). | Verificado contra o PAS com propofol e remifentanil juntos (diferença < 10⁻¹² mmHg). **Pendente:** conferir os 8 valores nas tabelas do artigo. |
| 14 | Variabilidade do remifentanil: ω² de V1, V2, V3, CL, Q2 e Q3 de Minto; efeito máximo na RPT (ω² 0,449) e inclinações na FC e no VS (covariância publicada) do Su. | PAS `pk_models.py` L375–L389; `pd_models.py` L1007, L1010, L1077–L1081. No PAS, os 3 do Su são **aditivos** (valor + η), e não × exp(η). | Testado (média e covariância). **Pendente:** conferir no artigo se são mesmo aditivos. |
| 15 | Com remifentanil, a farmacocinética do propofol usa os termos de Eleveld para opioide junto (V3 e depuração menores com a idade). | Eleveld 2018 (θ11 e θ13); PAS `pk_models.py` L310–L313. Por isso o BIS mínimo muda de 0 para qualquer taxa de remifentanil acima de 0. | Verificado contra o PAS. |
| 16 | O remifentanil **não entra no BIS**. | No modelo de BIS de Eleveld 2018, o termo do remifentanil é zero (PAS `pd_models.py` L224: `c50r = 0`). O PAS tem modelos de interação para o BIS (ex.: Bouillon 2004), mas trocar o modelo do BIS é outra decisão. | **Pendente:** decidir se usamos um modelo de interação para o BIS. |
| 17 | Entrada do remifentanil: **infusão contínua** (µg/kg/min), começando um tempo antes do propofol. | As bombas do VitalDB registram a taxa; o médico escolhe a taxa e o início. Taxa e início ficam fixos durante a simulação. | **Pendente:** opção de concentração-alvo (TCI), como as bombas do VitalDB usam. |
| 18 | Valores iniciais na plataforma: **0,12 µg/kg/min**, começando **9 min** antes do propofol. | Medianas das 21 induções do VitalDB do PBL_6: taxa média nos 20 min após o início do propofol 0,118 µg/kg/min (de 0,031 a 0,458); início do remifentanil 535 s antes do propofol (de 41 a 4182 s). Com esses valores, a Cp de Minto no início do propofol é 3,25 ng/ml no paciente de exemplo; nos 21 casos, a mediana foi 3,3 ng/ml. | — |
| 19 | PAM e FC basais: **antes de qualquer fármaco** (antes do remifentanil). | O modelo do Su parte do basal sem fármaco. No VitalDB do PBL_6, a PAM basal é medida 60 s antes do propofol, já com o remifentanil correndo. | **Pendente:** usar a PAM antes do remifentanil ao comparar com o VitalDB. |

Não usados (e por quê):
- **Variante "VitalDB" do Su no PAS** (`pd_models.py` L980–L987: `int_hr` −0,097, `fb` −0,5, efeito máximo do propofol na RPT −0,03). Foi ajustada junto com um modelo de estímulos (intubação e cirurgia) que não temos; sem ele, os valores não valem sozinhos.
- **Calibração da Beatriz no PBL_6** (`resultados/comparacao_propofol_remifentanil/`): é do modelo antigo (Windkessel, Schnider), não do Su.
- **`fator_fc_remifentanil` do PBL_6:** o próprio código diz que é uma hipótese exploratória, não validada.
- **Taxas acima de 0,5 µg/kg/min** no gráfico pressão × dose de remifentanil: a maior taxa média nos 21 casos foi 0,46.

Contagens feitas no VitalDB (6388 casos), para ajudar nessa decisão:
- 3512 casos com bomba de propofol; **3509** deles também com remifentanil. Numa amostra de 76 casos, o remifentanil começou sempre antes ou junto com o propofol (mediana: 30 s antes).
- 80 casos de anestesia inalatória com propofol em bolus manual, sem fentanil nem midazolam. O horário do bolus não é registrado (só a dose total) e só 6 têm pressão invasiva antes da indução.
