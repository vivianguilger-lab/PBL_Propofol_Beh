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
| 2 | Variação de Rp, VS e FC durante a anestesia (efeito do propofol e barorreflexo): estrutura e valores de Su 2023, como **complemento** do Windkessel (item 20). | Su H et al. *Br J Anaesth* 2023;131:222–33; implementação de referência no PAS (`pd_models.py` L890). Valores estimados em voluntários saudáveis: ponto de partida, não verdade final. | Verificado contra o PAS (diferença < 10⁻¹² mmHg). **Pendente:** conferir os 12 valores nas tabelas do artigo; o PAS declara `age_max_sv = 0.033` mas usa 0,0333; significado do "termo dependente do tempo". |
| 3 | ~~O VS basal não é usado.~~ **Substituída pelo item 20:** o VS basal agora é usado pelo Windkessel (VS0 = DC0/FC0). | No modelo de Su sozinho, a PAM não depende do VS basal (tudo em fração do basal). No Windkessel, o VS em mL define a onda de fluxo de cada batimento, a PAS e a PAD. | Resolvida. |
| 4 | Farmacocinética: Eleveld 2018, concentração arterial; **com os termos de opioide quando há remifentanil** (item 15). | Eleveld DJ et al. *Br J Anaesth* 2018;120:942–59. Mesmo grupo do Su; usado junto com o Su no PAS; população ampla; publica a variabilidade. | Verificado contra o PAS (diferença < 10⁻¹³ µg/ml). **Pendente:** qual PK o Su usou para calcular a Cp. |
| 5 | Entrada por concentração-alvo: alvo no **plasma**; taxa calculada com o paciente **típico**; sem limite de bomba por padrão. | O Su usa a Cp, e o alvo no sítio efetor exigiria um ke0. Uma bomba real calcula com o paciente típico. O limite da bomba depende do equipamento. | Testado (`taxa_alvo`). |
| 6 | Tempo em segundos em todo o código e nos gráficos. | Pedido do grupo, para facilitar a análise. | — |

## Windkessel (árvore arterial) — o centro do modelo

| # | Decisão | Por quê / fonte | Situação |
|---|---|---|---|
| 20 | A pressão arterial é calculada por um **Windkessel de 4 elementos** (forma paralela: Rc em paralelo com L, em série com C em paralelo com Rp), batimento a batimento. O modelo de Su entra como complemento: dá, a cada batimento, as frações do basal de Rp (RPT), VS e FC; a PAM do batimento volta ao Su como feedback do barorreflexo. | É o modelo da disciplina (Modelagem do Sistema Cardiovascular) e o modelo do projeto desde o início (PBL_6). Equações e parâmetros de Stergiopulos N et al. *Am J Physiol* 1999;276:H81–8; implementação de referência windkessel-id (Pigot, licença MIT), já portada no PBL_6. O Su elimina os valores ilustrativos de variação dos parâmetros que a versão anterior usava. | Testado: em regime, PAM do Windkessel = Rp × DC = PAM do Su; diferença máxima nos transitórios 0,43 mmHg (`runtests.jl`). |
| 21 | Discretização exata (a matriz do sistema é diagonal), passo de **5 ms**; o batimento começa no estado periódico x0 = xₙ⁰/(1 − aⁿ). | windkessel-id (`solvex0`); PBL_6 `src/windkessel.jl`. 5 ms é a amostragem do batimento de Stergiopulos. | Testado (o batimento basal se repete igual). |
| 22 | Onda de fluxo aórtico: o batimento tipo A de Stergiopulos (`dados/stergiopulos1999/aflow.csv`, 166 amostras), esticado para 60/FC e escalado para que o fluxo médio seja **exatamente** VS × FC/60. | Forma publicada e usada no windkessel-id. A escala exata garante a PAM basal igual à medida. | Testado. **Limitação:** a forma é a mesma para todos os pacientes. |
| 23 | Rc = 0,056 mmHg·s/mL e L = 0,0051 mmHg·s²/mL fixos, de Stergiopulos (1999). | Sem a onda de pressão invasiva, Rc e L não são identificáveis pelos dados do monitor. | **Pendente:** ajustar a partir da pressão invasiva do VitalDB, como no windkessel-id. |
| 24 | DC basal: o medido, se houver; se não, **3,0 L/min/m² × superfície corporal de Mosteller**. VS0 = DC0/FC0; Rp = PAM/DC0. | Índice cardíaco normal em repouso, como no PBL_6. Mosteller RD. *N Engl J Med* 1987;317:1098. Rp vem da lei de Ohm da circulação (no regime médio, C não passa fluxo e L curto-circuita Rc). | — |
| 25 | C identificada por bissecção (0,05 a 20 mL/mmHg) para que a **PAD** simulada seja a medida. | O Windkessel dá a pressão **central**. A PAM e a PAD quase não mudam da aorta ao braço; a PAS do manguito é amplificada na periferia (Agabiti-Rosei E et al. *Hypertension* 2007;50:154–60). Ajustar pela PAS do braço daria uma C errada. | Testado. Por isso a PAS central mostrada é menor que a do braço (20 a 31 mmHg nos três exemplos). |
| 26 | C constante durante a anestesia. | Su não descreve efeito do propofol na complacência; não há valor publicado para usar. | **Limitação.** |

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
- **Calibração da Beatriz no PBL_6** (`resultados/comparacao_propofol_remifentanil/`): calibra os fatores ilustrativos da versão anterior (sensibilidade e variação dos parâmetros do Windkessel com Schnider), que foram trocados pelos valores publicados do Su; não vale para as equações atuais.
- **`fator_fc_remifentanil` do PBL_6:** o próprio código diz que é uma hipótese exploratória, não validada.
- **Taxas acima de 0,5 µg/kg/min** no gráfico pressão × dose de remifentanil: a maior taxa média nos 21 casos foi 0,46.

Contagens feitas no VitalDB (6388 casos), para ajudar nessa decisão:
- 3512 casos com bomba de propofol; **3509** deles também com remifentanil. Numa amostra de 76 casos, o remifentanil começou sempre antes ou junto com o propofol (mediana: 30 s antes).
- 80 casos de anestesia inalatória com propofol em bolus manual, sem fentanil nem midazolam. O horário do bolus não é registrado (só a dose total) e só 6 têm pressão invasiva antes da indução.
