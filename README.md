# Antena Patch Retangular — Projeto Analítico (Balanis)

Script Python interativo para cálculo analítico completo de antenas patch retangulares com alimentação por microlinha e *inset feed*, baseado nas equações clássicas de Balanis (*Antenna Theory*, Cap. 14).

## O que o script calcula

A partir de frequência de operação, permissividade relativa do substrato, espessura do substrato e impedância característica desejada, o script resolve as 12 equações do projeto:

- Largura e comprimento do patch (`W`, `L`)
- Constante dielétrica efetiva (`εeff`)
- Comprimento efetivo e extensão por *fringing* (`Leff`, `ΔL`)
- Parâmetros da alimentação: largura da microlinha (`W0`), espaçamento do inset (`x0`)
- Condutância do slot e condutância mútua (`G1`, `G12`) via integração numérica
- Resistência de entrada (`Rin`) e profundidade/comprimento do inset (`y0`, `L0`)

A largura da microlinha (`W0`) é obtida resolvendo numericamente a equação de impedância característica (coeficiente 0,667, conforme Balanis — não 0,337 como aparece com erro tipográfico em alguns artigos de referência).

## Uso

```bash
python projeto_antena_patch.py
```

O programa pede interativamente:
- Frequencia f [GHz]:
- Constante dieletrica er:
- Espessura h [mm]:
- Impedancia alvo Zc [Ohm] [50]:

Pressionar Enter no campo de impedância usa o padrão de 50 Ω.

## Saída

- Tabela formatada no terminal com todos os parâmetros calculados
- Arquivo `resultados_antena.csv`, com cada execução adicionada como nova linha (histórico cumulativo de projetos)

## Dependências

```bash
pip install numpy scipy
```

## Referência

Balanis, C. A. *Antenna Theory: Analysis and Design*. Cap. 14 — Microstrip Antennas.
