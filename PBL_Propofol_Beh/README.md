# PropofolWK

Previsão da pressão arterial durante a indução anestésica com **propofol**
(e, se o médico escolher, **remifentanil**),
apresentada como **faixas** — e não como uma linha única — com a indicação
dos períodos em que o paciente permanece em cada nível de risco de hipotensão.

## Princípios

- Nenhum valor fixo de artigo é tratado como verdade para todos os pacientes:
  os valores vêm de pacientes reais (VitalDB) e são mostrados como faixas.
- Toda estrutura de modelo usada tem base publicada, e toda decisão de método
  fica registrada em `docs/`, com a fonte de cada valor.
- O que ainda não tem base fica fora do modelo até ser validado.

## Entrada (escolha do médico)

Propofol:
1. **Dose e velocidade de infusão** (mg/kg e mg/kg/h); ou
2. **Concentração-alvo** (µg/ml).

Remifentanil (opcional): **infusão contínua** (µg/kg/min) e quantos minutos antes
do propofol ela começou.

Dados do paciente: idade, sexo, peso, altura, ASA e sinais pré-indução
(PAS, PAD, PAM, FC).

## Saída

- Faixa da PAM prevista ao longo do tempo (em segundos após o início do bolus),
  desenhada sobre as zonas de risco.
- Para cada intervalo de tempo: entre quais valores a PAM fica e em que nível de
  risco o paciente permanece.

## Bases

- **Dados:** VitalDB — Lee HC et al. *Sci Data* 2022;9:279.
- **Efeito do propofol na circulação (estrutura):** Su H et al. *Br J Anaesth*
  2023;131:222–33; implementação de referência:
  [Python Anesthesia Simulator](https://github.com/AnesthesiaSimulation/Python_Anesthesia_Simulator).
- **Farmacocinética do remifentanil:** Minto CF et al. *Anesthesiology*
  1997;86:10–23 (o mesmo modelo das bombas do VitalDB; ver `docs/decisoes.md`).
- **Circulação:** Windkessel de 4 elementos — Stergiopulos N et al.
  *Am J Physiol* 1999;276:H81–8; implementação de referência:
  [windkessel-id](https://github.com/hpigot/windkessel-id).
- **Seleção de casos do VitalDB (referência):**
  [PAS_vs_vitalDB](https://github.com/AnesthesiaSimulation/PAS_vs_vitalDB).

## Situação

Em construção...
