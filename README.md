# Analisis-Numerico---TP2
Métodos numéricos aplicados a la Ingeniería de procesos

Uno de los procesos más comunes de la industria metalúrgica es el tratamiento térmico de
aceros. Las múltiples aplicaciones de esta aleación Fe-C (automotriz, petróleo, tuberías,
perfiles, etc.) requieren propiedades mecánicas muy diversas.
En este sentido, los tratamientos térmicos suelen ser una herramienta versátil y rentable para
lograr un objetivo de resistencia mecánica y/o tenacidad.
Para una misma composición química es posible modificar las propiedades mecánicas del acero
sometiendo el material a una serie de calentamientos y enfriamientos sucesivos.
Uno de los tratamientos térmicos de gran aplicación en la industria es el Temple y Revenido,
que consiste en las siguientes etapas:
1) Austenizado: Calentamiento hasta una temperatura de aproximadamente 900°C. En este
proceso el material modifica completamente su microestructura (transformación alotrópica) y
cambia a la fase Austenita.
2) Temple: Enfriamiento brusco hasta temperatura ambiente, generalmente utilizando agua. Se
obtiene una nueva fase llamada Martensita. En esta condición el material tiene alta resistencia
mecánica pero a la vez es muy frágil (baja tenacidad).
3) Revenido: Nuevo calentamiento hasta una temperatura objetivo. El rango habitual es 500-
700°C. A medida que el material se calienta la Martensita disminuye su resistencia mecánica
ganando tenacidad. La temperatura final definirá la combinación final de propiedades
mecánicas.

El presente trabajo se focalizará en modelar la evolución temporal de la tercera etapa dicho
tratamiento térmico. El material a revenir serán tubos de acero destinados a la industria
petrolera. La norma API es la encargada de definir los requerimientos de propiedades
mecánicas según la profundidad y las condiciones particulares del pozo. En este contexto es
necesario contar con un proceso que garantice el producto final dentro de rango de propiedades
mecánicas.
Para lograr un proceso continuo de calentamiento, los tubos avanzan dentro de un horno a
velocidad constante. Cada cierto intervalo de tiempo (cadencia) un tubo sale del horno a la
temperatura objetivo, a la vez un nuevo tubo entra a temperatura ambiente.
Por otra parte, el horno debe tener la suficiente versatilidad para calentar tubos de distintas
dimensiones, así como también lograr diferentes objetivos de temperatura. Por este motivo los
hornos son diseñados con varias zonas, pudiendo establecer temperaturas independientes en
cada una.
