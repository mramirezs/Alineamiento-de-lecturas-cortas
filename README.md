# Alineamiento-de-lecturas-cortas

## **Objetivo de la Sesión**
- Aprender el concepto y la práctica del alineamiento de lecturas cortas.
- Usar herramientas comunes como **BWA** y **Bowtie2** para trabajar con datos reales de *Ebola*.
- Comprender los formatos estándar como **SAM/BAM** y analizar los resultados.
  
---

## **1. Introducción Teórica**

### **1.1 ¿Qué es el alineamiento de lecturas cortas?**
- **Definición:** Es el proceso de comparar lecturas de ADN o ARN (obtenidas mediante tecnologías de secuenciación como Illumina) con un genoma de referencia (una secuencia de ADN representativa que se utiliza como estándar para comparar otras secuencias genómicas).
- **Importancia:** Permite identificar variantes genéticas, analizar la expresión génica, entre otros estudios genómicos.

### **1.2 Conceptos clave**
1. **Reads (Lecturas):** Fragmentos cortos de ADN/ARN (150-300 bp) generados por secuenciación.
2. **Genoma de referencia:** Secuencia completa de un organismo utilizada como punto de comparación.
3. **Indexación:** Preprocesamiento del genoma de referencia para acelerar el proceso de búsqueda de coincidencias.
4. **SAM/BAM:** Formatos para almacenar resultados de alineación:
   - **SAM:** Texto plano.
   - **BAM:** Formato binario comprimido.

### **1.3 Herramientas a utilizar**
1. **BWA (Burrows-Wheeler Aligner):** Optimizado para datos de Illumina, soporta lecturas de longitud media.
2. **Bowtie2:** Rápido y flexible, útil para alineaciones con diferencias significativas.
3. **Samtools:** Procesamiento y análisis de resultados en formatos SAM/BAM.

---

## **2. Actividades Prácticas**

### **2.1 Preparación del entorno bioinfomático para el proceso de datos**
1. **Instalar herramientas necesarias:**
   ```bash
   sudo apt install bwa bowtie2 samtools
   pip install bio --upgrade
   ```
2. **Descargar de los datos a procesar y analizar:**
   - Genoma de referencia (*Ebola*, cepa de 1976):
     ```bash
     mkdir -p refs
     bio fetch AF086833 --format fasta > refs/ebola_ref.fa
     ```
   - Datos de secuenciación (*SRR1972739* (Zaire ebolavirus genome, 2014), 10,000 lecturas):
     ```bash
     mkdir -p data
     fastq-dump -X 10000 --split-files SRR1972739 -O data
     ```

---

### **2.2 Creación del Índice del Genoma**

1. **¿Qué es un índice?**
   - Un índice es una estructura de datos que permite buscar coincidencias en el genoma de referencia de manera eficiente.

2. **Crear índice con BWA:**

   ```bash
   bwa index refs/ebola_ref.fa
   ls refs/
   ebola_ref.fa  ebola_ref.fa.amb  ebola_ref.fa.ann  ebola_ref.fa.bwt  ebola_ref.fa.pac  ebola_ref.fa.sa
   ```
   
4. **Crear índice con Bowtie2:**

   ```bash
   bowtie2-build refs/ebola_ref.fa refs/ebola_ref
   ebola_ref.1.bt2  ebola_ref.3.bt2  ebola_ref.fa ebola_ref.rev.2.bt2
   ebola_ref.2.bt2  ebola_ref.4.bt2  ebola_ref.rev.1.bt2
   ```
   
5. **Verificar los archivos generados:**
   ```bash
   ls refs/
   ebola_ref.1.bt2  ebola_ref.3.bt2  ebola_ref.fa      ebola_ref.fa.ann  ebola_ref.fa.pac  ebola_ref.rev.1.bt2
   ebola_ref.2.bt2  ebola_ref.4.bt2  ebola_ref.fa.amb  ebola_ref.fa.bwt  ebola_ref.fa.sa   ebola_ref.rev.2.bt2
   ```

---

### **2.3 Alineación de Lecturas**
#### **1. Alineación con BWA**
- **Modo paired-end:** Alinea las dos lecturas (*forward* y *reverse*) simultáneamente:
  ```bash
  bwa mem refs/ebola_ref.fa data/SRR1972739_1.fastq data/SRR1972739_2.fastq > bwa_output.sam
  ```
- **Verificar resultados iniciales:**
  ```bash
  head -n 20 bwa_output.sam
  ```

#### **2. Alineación con Bowtie2**
- **Modo paired-end:**
  ```bash
  bowtie2 -x refs/ebola_ref -1 data/SRR1972739_1.fastq -2 data/SRR1972739_2.fastq -S bowtie2_output.sam
  10000 reads; of these:
  10000 (100.00%) were paired; of these:
    5086 (50.86%) aligned concordantly 0 times
    4914 (49.14%) aligned concordantly exactly 1 time
    0 (0.00%) aligned concordantly >1 times
    ----
    5086 pairs aligned concordantly 0 times; of these:
      827 (16.26%) aligned discordantly 1 time
    ----
    4259 pairs aligned 0 times concordantly or discordantly; of these:
      8518 mates make up the pairs; of these:
        7463 (87.61%) aligned 0 times
        1055 (12.39%) aligned exactly 1 time
        0 (0.00%) aligned >1 times
  62.69% overall alignment rate
  ```
  
- **Verificar resultados iniciales:**
  ```bash
  head -n 20 bowtie2_output.sam
  ```

---

### **2.4 Análisis de Resultados**
#### **1. Comprender el formato SAM**
- Cada línea contiene:
  - Identificador de la lectura.
  - Coordenadas del alineamiento.
  - Calidad del mapeo.
    
- Ejemplo de línea:
  ```
  SRR1972739.1   83   gi|10141003|gb|AF086833.2|   15684   60   69M32S   =   15600   -153   TTTAGATTTAAC
  ```

#### **2. Estadísticas básicas con Samtools**
- Generar estadísticas:
  ```bash
  samtools flagstat bwa_output.sam
  samtools flagstat bowtie2_output.sam
  ```

#### **3. Conversión a BAM y visualización**
- Convertir a formato BAM:
  
  ```bash
  samtools view -b bwa_output.sam > bwa_output.bam
  samtools view -b bowtie2_output.sam > bowtie2_output.bam
  ```
- Visualizar alineaciones:
  ```bash
  samtools tview bwa_output.bam refs/ebola_ref.fa
  ```

---

## **3. Conclusión y Tareas**

### Resumen de la sesión:
- Aprendimos los pasos básicos del alineamiento de lecturas cortas.
- Usamos BWA y Bowtie2 para realizar alineaciones.
- Analizamos y visualizamos los resultados en formatos SAM/BAM.

### Tareas sugeridas:
1. Experimentar con datos propios o nuevos:
   - Descargar otro genoma de referencia.
   - Realizar alineaciones con lecturas más largas o con mayor número de datos.
     
2. Explorar parámetros avanzados:
   ```bash
   bwa mem -t 4 refs/ebola_ref.fa data/SRR1972739_1.fastq > advanced_output.sam
   bowtie2 --very-sensitive -x refs/ebola_ref -1 data/SRR1972739_1.fastq > advanced_output.sam
   ```
3. Investigar más sobre los formatos SAM/BAM y su manejo con Samtools.
