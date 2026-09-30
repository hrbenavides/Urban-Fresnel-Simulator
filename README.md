# Urban Fresnel Simulator

Simulador tridimensional de planificación, análisis de factibilidad y cálculo de presupuesto de potencia (*Link Budget*) para radioenlaces urbanos, basado en el análisis del **elipsoide de Fresnel** y la **difracción por obstáculos** según la recomendación **UIT-R P.526**.

Herramienta de ingeniería diseñada para escenarios metropolitanos densos, donde los edificios y las estructuras civiles representan el mayor desafío de propagación.

**Autor:** Univ. Henry Rafael Benavides Gutierrez  
**Institución:** Universidad Mayor de San Andrés — Facultad de Ingeniería

---

## 📥 Descarga

¿Solo quieres usar el programa? **Descarga el ejecutable desde la sección [Releases](../../releases)**.

- ✅ No requiere instalar Python.
- ✅ No requiere instalar dependencias.
- ✅ Doble clic y a trabajar.
- 🖥️ Compatible con Windows 10 / 11 (64-bit).

---

## ✨ Características principales

- **Análisis de superficie de techo adaptativo:** escaneo tridimensional de tapas y aristas de las estructuras para hallar el punto crítico de aproximación macroscópica (*t\**).
- **Motor geométrico 3D real:** soporta rotación axial independiente de las estructuras y posicionamiento de antenas en cualquiera de las 4 esquinas del tejado o su eje central.
- **Cálculo de potencia integrado:** evalúa dinámicamente la Pérdida por Propagación en el Espacio Libre (FSPL), acoplada a castigos por difracción de la UIT-R ante la invasión del elipsoide.
- **Análisis multi-obstáculo:** soporte para hasta **6 obstáculos urbanos** configurables en posición, altura, ancho, offset lateral y rotación.
- **Presupuesto de potencia completo:** potencia Tx (dBm / W / kW), ganancias de antena (dB / dBi / dBd), sensibilidad del receptor (dBm / W / kW / µW / nW) y margen de seguridad configurable.
- **Interfaz estilo CAD:** renderizado tridimensional interactivo con estética oscura (*Blender Dark*), modo Rayos X (estructuras translúcidas) y cotas dinámicas en el lienzo.
- **Veredicto automático de factibilidad:** evaluación combinada de geometría (despeje al 60 %) y potencia (margen respecto a la sensibilidad).

---

## 📊 Fundamentos matemáticos

**Radio de la *n*-ésima zona de Fresnel:**

$$r_n = \sqrt{\frac{n \cdot \lambda \cdot d_1 \cdot d_2}{D_{\text{real}}}}$$

**Atenuación por difracción de filo de cuchillo (UIT-R P.526):**

$$\nu = H \sqrt{\frac{2 (d_1 + d_2)}{\lambda \, d_1 \, d_2}}$$

$$A_{\text{dif}} = 6.9 + 20 \log_{10} \left( \sqrt{(\nu - 0.1)^2 + 1} + \nu - 0.1 \right)$$

**Presupuesto de potencia:**

$$P_{rx}\ (\text{dBm}) = P_{tx} + G_{tx} + G_{rx} - FSPL - A_{\text{dif}} - M_{\text{seg}}$$

donde $FSPL = 20\log_{10}(d) + 20\log_{10}(f) - 147.55$.

---

## 🖥️ Requisitos

- **Windows 10 / 11 (64-bit)** para el ejecutable.
- **Python 3.10 o superior** (probado en 3.14) si se ejecuta desde el código fuente.
- **Tkinter** (incluido en la instalación estándar de Python en Windows).

Dependencias de Python:

- `numpy`
- `matplotlib`

---

## 🚀 Ejecutar desde el código fuente

Si quieres modificar el código o contribuir:

```bash
# 1. Clonar el repositorio
git clone https://github.com/hrbenavides/Urban-Fresnel-Simulator.git
cd Urban-Fresnel-Simulator

# 2. Crear y activar un entorno virtual
python -m venv venv
venv\Scripts\activate

# 3. Instalar dependencias
pip install -r requirements.txt

# 4. Ejecutar
python fresnel15.py