"""
Programa interativo: informe f, er, h e Zc no terminal.
Equacoes 1-12 conforme artigo + Balanis.
Saida: tabela no terminal + arquivo CSV.
"""

import csv
import os
import numpy as np
from scipy import integrate
from scipy.special import j0
from scipy.optimize import brentq


def ler_float(msg, padrao=None):
    """Le um float do terminal. Aceita virgula ou ponto. Enter usa o padrao."""
    while True:
        txt = input(msg).strip().replace(",", ".")
        if txt == "" and padrao is not None:
            return padrao
        try:
            return float(txt)
        except ValueError:
            print("  >> Valor invalido. Tente novamente.")


def projetar(f, er, h, Zc):
    """Calcula todos os parametros da antena patch. Retorna dict."""
    c = 3e8

    # Eq. 1 — Largura do patch
    W = (c / (2 * f)) * np.sqrt(2 / (er + 1))

    # Eq. 2 — Constante dieletrica efetiva
    eeff = (er + 1) / 2 + (er - 1) / 2 * (1 / np.sqrt(1 + 12 * h / W))

    # Eq. 4 — Comprimento efetivo
    Leff = c / (2 * f * np.sqrt(eeff))

    # Eq. 5 — Extensao de comprimento dL (fringing)
    dL = 0.412 * h * ((eeff + 0.3) * (W / h + 0.264)) / \
         ((eeff - 0.258) * (W / h + 0.8))

    # Eq. 3 — Comprimento real do patch
    L = Leff - 2 * dL

    # Alimentacao
    lo = c / f                  # comprimento de onda no ar
    k0 = 2 * np.pi / lo         # numero de onda

    # Eq. 6 — Espacamento antena-linha
    x0 = lo / 100

    # Eq. 7 — Largura da microlinha para Zc (Balanis, coef. 0.667)
    def eq7_residual(W0_val):
        w0sobH = W0_val / h
        return (120 * np.pi / np.sqrt(eeff)) / \
               (w0sobH + 1.393 + 0.667 * np.log(w0sobH + 1.444)) - Zc
    W0 = brentq(eq7_residual, 0.01e-3, 50e-3)

    # Eq. 10 — Condutancia de um slot
    G1 = (1 / 120) * (W / lo)

    # Eq. 11 — Condutancia mutua G12 (integracao numerica)
    def integrand_G12(theta):
        cos_t = np.cos(theta)
        sin_t = np.sin(theta)
        arg = k0 * W / 2 * cos_t
        frac = np.where(np.abs(cos_t) < 1e-12,
                        k0 * W / 2,
                        np.sin(arg) / cos_t)
        return (frac ** 2) * j0(k0 * L * sin_t) * sin_t ** 3
    G12_integral, _ = integrate.quad(integrand_G12, 0, np.pi)
    G12 = G12_integral / (120 * np.pi ** 2)

    # Eq. 9 — Resistencia de entrada em y=0 (sinal + : distribuicao par)
    Rin_0 = 1 / (2 * (G1 + G12))

    # Eq. 8 — Posicao do inset y0
    y0 = (L / np.pi) * np.arccos(np.sqrt(Zc / Rin_0))

    # Eq. 12 — Comprimento da microlinha
    L0 = 2 * y0

    return {
        "f_GHz": f / 1e9,
        "er": er,
        "h_mm": h * 1e3,
        "Zc_ohm": Zc,
        "W_mm": W * 1e3,
        "eeff": eeff,
        "Leff_mm": Leff * 1e3,
        "dL_mm": dL * 1e3,
        "L_mm": L * 1e3,
        "lambda0_mm": lo * 1e3,
        "x0_mm": x0 * 1e3,
        "W0_mm": W0 * 1e3,
        "G1_S": G1,
        "G12_S": G12,
        "Rin0_ohm": Rin_0,
        "y0_mm": y0 * 1e3,
        "L0_mm": L0 * 1e3,
    }


def imprimir(r):
    print("\n" + "=" * 60)
    print("  ANTENA PATCH — RESULTADOS")
    print("=" * 60)
    print(f"\n  [ENTRADA]  f = {r['f_GHz']} GHz | er = {r['er']} | "
          f"h = {r['h_mm']} mm | Zc = {r['Zc_ohm']} Ohm")
    print(f"\n  [PATCH]")
    print(f"  W    = {r['W_mm']:.4f} mm       (largura da patch)")
    print(f"  eeff = {r['eeff']:.4f}           (const. dieletrica efetiva)")
    print(f"  Leff = {r['Leff_mm']:.4f} mm       (comprimento efetivo)")
    print(f"  dL   = {r['dL_mm']:.4f} mm        (correcao de comprimento)")
    print(f"  L    = {r['L_mm']:.4f} mm       (comprimento da patch)")
    print(f"\n  [ALIMENTACAO]")
    print(f"  lambda_0 = {r['lambda0_mm']:.2f} mm        (comprimento de onda)")
    print(f"  x0       = {r['x0_mm']:.4f} mm        (espacamento)")
    print(f"  W0       = {r['W0_mm']:.4f} mm        (largura da microlinha)")
    print(f"  G1       = {r['G1_S']:.6e} S    (condutancia do slot)")
    print(f"  G12      = {r['G12_S']:.6e} S    (condutancia mutua)")
    print(f"  Rin(y=0) = {r['Rin0_ohm']:.4f} Ohm      (resistencia de entrada)")
    print(f"  y0       = {r['y0_mm']:.4f} mm        (profundidade de insercao)")
    print(f"  L0       = {r['L0_mm']:.4f} mm        (comprimento da microlinha)")
    print("=" * 60)


def salvar_csv(r, nome="resultados_antena.csv"):
    """Anexa o resultado ao CSV (cria com cabecalho se nao existir)."""
    existe = os.path.isfile(nome)
    with open(nome, "a", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=list(r.keys()))
        if not existe:
            writer.writeheader()
        writer.writerow(r)
    print(f"\n  >> Resultados salvos em: {nome}")


def main():
    print("=" * 60)
    print("  PROJETO DE ANTENA PATCH RETANGULAR — IC 2026")
    print("  (Enter aceita o valor padrao entre colchetes)")
    print("=" * 60)

    f  = ler_float("\n  Frequencia f [GHz]: ") * 1e9
    er = ler_float("  Constante dieletrica er : ")
    h  = ler_float("  Espessura h [mm]: ") * 1e-3
    Zc = ler_float("  Impedancia alvo Zc [Ohm] [50]: ", padrao=50.0)

    r = projetar(f, er, h, Zc)
    imprimir(r)
    salvar_csv(r)


if __name__ == "__main__":
    main()