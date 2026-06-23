# Multi-Model Fresnel Zone & RF Link Budget Simulator

Este repositorio contiene una suite avanzada de simulación tridimensional y bidimensional para la planificación, análisis de factibilidad y cálculo de presupuesto de potencia (*Link Budget*) en radioenlaces de telecomunicaciones. El software ha sido desarrollado con un enfoque de alta fidelidad visual y precisión matemática milimétrica, sirviendo como herramienta de ingeniería tanto para entornos urbanos complejos como para enlaces rurales de larga distancia.

**Autor:** Univ. Henry Rafael Benavides Gutierrez  
**Institución:** Universidad Mayor de San Andrés — Facultad de Ingeniería  

---

## 🚀 Características Principales

La suite se divide en dos módulos core especializados e independientes:

### 1. Módulo Urbano 3D (`fresnel15.py`)
Diseñado para escenarios metropolitanos densos donde las estructuras civiles y la orientación de los edificios representan el mayor desafío de propagación.
* **Análisis de Superficie de Techo Adaptativo:** Escanea mediante matrices tridimensionales las tapas y aristas de las estructuras para hallar el punto crítico de aproximación macroscópica ($t^*$).
* **Motor Geométrico Tridimensional Real:** Soporta la rotación axial independiente de las estructuras y el posicionamiento de antenas en cualquiera de las 4 esquinas del tejado o su eje central.
* **Cálculo de Potencia Integrado:** Evalúa de forma dinámica la Pérdida por Propagación en el Espacio Libre (FSPL) acoplada a castigos por difracción de la UIT-R ante la invasión del elipsoide.
* **Interfaz Estilo CAD:** Renderizado tridimensional interactivo con estética oscura (*Blender Dark*), soporte para modo Rayos X (estructuras transparentes) y cotas dinámicas en el lienzo.

---

## 📊 Fundamentos Matemáticos Utilizados

* **Radio de la $n$-ésima Zona de Fresnel:**
  $$r_n = \sqrt{\frac{n \cdot \lambda \cdot d_1 \cdot d_2}{D_{\text{real}}}}$$
* **Flecha de Curvatura Terrestre con Refracción Atmosférica:**
  $$h_{crv} = \frac{d_1 \cdot (D_{\text{total}} - d_1)}{2 \cdot K \cdot R_e}$$
* **Presupuesto de Potencia (Friis & UIT-R P.526):**
  $$P_{rx}\text{ (dBm)} = P_{tx} + G_{tx} + G_{rx} - FSPL - A_{dif}$$

---

## 🛠️ Requisitos e Instalación

El sistema está desarrollado íntegramente en Python 3 utilizando la biblioteca nativa `tkinter` para la interfaz gráfica, asegurando un entorno liviano y portable.

1. Clona este repositorio:
   ```bash
   git clone [https://github.com/hrbenavides/fresnel-link-simulator.git](https://github.com/hrbenavides/fresnel-link-simulator.git)
   cd fresnel-link-simulator
