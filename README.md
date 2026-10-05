<div align="center">

# 📡 Urban Fresnel Simulator

**Simulador profesional de enlaces de radiofrecuencia urbanos**  
Cálculo geométrico de zona de Fresnel, análisis multi-obstáculo y veredicto automático de viabilidad.

[![Descargar ejecutable](https://img.shields.io/badge/⬇️_Descargar-Ejecutable-blue?style=for-the-badge)](https://github.com/hrbenavides/Urban-Fresnel-Simulator/releases)


[![Python](https://img.shields.io/badge/Python-3.9+-3776AB?style=flat-square&logo=python&logoColor=white)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-green?style=flat-square)](LICENSE)
[![Estado](https://img.shields.io/badge/Estado-Activo-brightgreen?style=flat-square)]()

</div>

---

## 🌐 Descripción general

**Urban Fresnel Simulator** es una herramienta de escritorio desarrollada en Python para el **análisis y diseño de enlaces de telecomunicaciones en entornos urbanos**. Permite evaluar la viabilidad de un enlace punto a punto considerando obstáculos reales (edificios, terreno, vegetación) y aplicando los modelos de propagación recomendados por la UIT.

El simulador integra un motor geométrico propio, una interfaz gráfica tipo CAD y un sistema de análisis automático que entrega un **veredicto claro sobre la calidad del enlace**.

---

## ✨ Características principales

- 🧮 **Motor geométrico propio** — cálculo de la primera zona de Fresnel y su porcentaje de obstrucción.
- 🏙️ **Análisis multi-obstáculo** — soporta múltiples edificios y perfiles de terreno entre el transmisor y el receptor.
- 🖥️ **Interfaz gráfica tipo CAD** — dibujo interactivo de obstáculos sobre el perfil del enlace.
- 📊 **Modelo de difracción UIT-R P.526** — atenuación por filo de cuchillo aplicada automáticamente.
- ⚡ **Presupuesto de potencia** — cálculo de FSPL, potencia recibida y margen del enlace.
- ✅ **Veredicto automático** — clasifica el enlace como *viable*, *marginal* o *no viable*.
- 🎨 **Modo Rayos X** — visualización esquemática del trazado del enlace.
- 💾 **Exportación de resultados** — guarda reportes en formato de texto e imagen.

---

## 🛠️ Tecnologías utilizadas

| Capa | Tecnología |
|------|------------|
| Lenguaje | Python 3.9+ |
| Interfaz gráfica | PyQt5 / PySide2 |
| Cálculo numérico | NumPy |
| Gráficos | Matplotlib |
| Empaquetado | PyInstaller |
| Distribución | GitHub Releases |

---

## 🚀 Instalación y uso

### Opción 1 — Descargar el ejecutable (recomendado)

1. Ve a la sección de **[Releases](https://github.com/hrbenavides/Urban-Fresnel-Simulator/releases)**.
2. Descarga el archivo correspondiente a tu sistema operativo.
3. Descomprime y ejecuta `UrbanFresnel.exe` (Windows) o el binario equivalente.
4. No requiere instalación de Python ni dependencias.

### Opción 2 — Ejecutar desde el código fuente

```bash
# 1. Clonar el repositorio
git clone https://github.com/hrbenavides/Urban-Fresnel-Simulator.git
cd Urban-Fresnel-Simulator

# 2. Crear y activar un entorno virtual
python -m venv venv

# Windows:
venv\Scripts\activate

# Linux / macOS:
source venv/bin/activate

# 3. Instalar dependencias
pip install -r requirements.txt

# 4. Ejecutar el simulador
python fresnel15.py
```

---

## 📐 Fundamentos matemáticos

El simulador implementa los modelos de propagación de la **Unión Internacional de Telecomunicaciones (UIT)** para el análisis de enlaces en entornos urbanos.

### Radio de la n-ésima zona de Fresnel

$$r_n = \sqrt{\frac{n \lambda d_1 d_2}{d_1 + d_2}}$$

donde:

- $r_n$ — radio de la n-ésima zona de Fresnel (m)
- $\lambda$ — longitud de onda (m)
- $d_1$, $d_2$ — distancias del obstáculo al transmisor y receptor (m)

### Atenuación por difracción — filo de cuchillo (UIT-R P.526)

$$A_d = 6.9 + 20 \log_{10}\left(\sqrt{(v - 0.1)^2 + 1} + v - 0.1\right)$$

donde $v$ es el parámetro de Fresnel-Kirchhoff:

$$v = h \sqrt{\frac{2 (d_1 + d_2)}{\lambda d_1 d_2}}$$

### Presupuesto de potencia del enlace

$$P_{rx} = P_{tx} + G_{tx} + G_{rx} - FSPL - A_d$$

donde el **espacio libre (FSPL)** se calcula como:

$$FSPL = 20 \log_{10}(d) + 20 \log_{10}(f) - 147.55$$

con $d$ en metros y $f$ en Hz.

---

## 📊 Interpretación del veredicto

El simulador clasifica el enlace según el porcentaje de obstrucción de la primera zona de Fresnel:

| Obstrucción | Veredicto | Interpretación |
|:---:|:---:|:---|
| **< 20 %** | ✅ **Viable** | Enlace con margen suficiente para operar. |
| **20 % – 40 %** | ⚠️ **Marginal** | Enlace operable pero con degradación notable. |
| **> 40 %** | ❌ **No viable** | Obstrucción crítica, requiere rediseño. |

---

## 📁 Estructura del proyecto

```
Urban-Fresnel-Simulator/
├── fresnel15.py           # Script principal del simulador
├── requirements.txt       # Dependencias del proyecto
├── LICENSE                # Licencia MIT
├── README.md              # Este archivo
├── docs/                  # Capturas y documentación adicional
│   └── screenshot-*.png
└── dist/                  # Ejecutables generados con PyInstaller
```

---

## 🤝 Contribuciones

Las contribuciones son bienvenidas. Si deseas mejorar el simulador:

1. Haz un **fork** del repositorio.
2. Crea una rama para tu mejora: `git checkout -b feature/nueva-funcionalidad`.
3. Realiza tus cambios y haz commit: `git commit -m "Añade nueva funcionalidad"`.
4. Haz push a tu rama: `git push origin feature/nueva-funcionalidad`.
5. Abre un **Pull Request**.

---

## 📝 Licencia

Este proyecto está bajo la licencia **MIT**. Consulta el archivo [LICENSE](LICENSE) para más detalles.

---

## 👨‍💻 Autor

**Henry Benavides**

- GitHub: [@hrbenavides](https://github.com/hrbenavides)
- Proyecto desarrollado como parte del curso de **Ingeniería Electrónica — UMSA**

---

<div align="center">

## 📥 ¿Listo para probarlo?

👉 **[Descargar la última versión aquí](https://github.com/hrbenavides/Urban-Fresnel-Simulator/releases)**

<br>

Hecho en Bolivia ❤️💛💚

</div>
