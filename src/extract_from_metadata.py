#!/usr/bin/env python3
import csv
import sys
import os

def main():
    csv_path = 'data/GSE213001/SraRunTable.csv'
        
    # Lista de las 12 muestras NDC (sanas) y 12 muestras IPF seleccionadas
    selected_gsms = {
        # NDC (Sanas)
        "GSM6568449", "GSM6568325", "GSM6568417", "GSM6568367",
        "GSM6568382", "GSM6568349", "GSM6568361", "GSM6568399",
        "GSM6568445", "GSM6568389", "GSM6568359", "GSM6568353",
        # IPF
        "GSM6568380", "GSM6568437", "GSM6568330", "GSM6568373",
        "GSM6568438", "GSM6568393", "GSM6568317", "GSM6568443",
        "GSM6568356", "GSM6568342", "GSM6568348", "GSM6568387"
    }
    
    results = []
    seen_gsm = set()
    
    with open(csv_path, mode='r', encoding='utf-8') as f:
        # skipinitialspace=True maneja comillas y espacios iniciales
        reader = csv.DictReader(f, skipinitialspace=True)
        # Limpiar espacios en los nombres de las columnas
        reader.fieldnames = [name.strip() for name in reader.fieldnames]
        
        for row in reader:
            gsm = row['Sample Name'].strip()
            if gsm in selected_gsms:
                if gsm not in seen_gsm:
                    seen_gsm.add(gsm)
                    results.append({
                        'diagnosis': row['diagnosis'].strip(),
                        'sample_name': gsm,
                        'donor_id': row['donorid'].strip(),
                        'age': row['AGE'].strip(),
                        'gender': row['gender'].strip(),
                        'lunglocation': row['lunglocation'].strip(),
                        'smokingstatus': row['smokingstatus'].strip(),
                        'severity': row['Severity'].strip(),
                        'run': row['Run'].strip()
                    })
                    
    # Ordenar los resultados por diagnóstico (NDC primero, luego IPF) y por nombre de muestra
    results.sort(key=lambda x: (x['diagnosis'], x['sample_name']))
    
    # Imprimir cabecera de TSV
    headers_list = ["Diagnosis", "Sample_Name", "Donor_ID", "Age", "Gender", "Location", "Smoking", "Severity", "Run"]
    print("\t".join(headers_list))
    
    # Imprimir filas
    for r in results:
        row_values = [
            r['diagnosis'],
            r['sample_name'],
            r['donor_id'],
            r['age'],
            r['gender'],
            r['lunglocation'],
            r['smokingstatus'],
            r['severity'],
            r['run']
        ]
        print("\t".join(row_values))

if __name__ == '__main__':
    main()
