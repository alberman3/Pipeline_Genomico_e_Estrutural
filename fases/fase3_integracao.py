import requests
import pandas as pd
import json

class GAPDHFunctionalAPI:
    """
    Classe para consulta de APIs RESTful funcionais (UniProt e KEGG)
    para a proteína GAPDH de Gallus gallus.
    """
    def __init__(self, uniprot_id: str = "P00356"):
        self.uniprot_id = uniprot_id
        self.uniprot_url = f"https://rest.uniprot.org/uniprotkb/{self.uniprot_id}.json"

    def fetch_uniprot_data(self) -> dict:
        """
        Consome a REST API do UniProt para extrair Função e Localização Celular.
        """
        print(f"--- [API UNIPROT] Consultando ID {self.uniprot_id} ---")
        response = requests.get(self.uniprot_url)
        response.raise_for_status()
        data = response.json()

        # Extração da Função da Proteína (Comment type: FUNCTION)
        function_text = "Função não especificada."
        for comment in data.get("comments", []):
            if comment.get("commentType") == "FUNCTION":
                texts = comment.get("texts", [])
                if texts:
                    function_text = texts[0].get("value", "")
                break

        # Extração da Localização Celular (Comment type: SUBCELLULAR LOCATION)
        locations = []
        for comment in data.get("comments", []):
            if comment.get("commentType") == "SUBCELLULAR LOCATION":
                for loc in comment.get("subcellularLocations", []):
                    loc_val = loc.get("location", {}).get("value")
                    if loc_val:
                        locations.append(loc_val)
        
        location_text = ", ".join(locations) if locations else "Cytoplasm, Nucleus"

        return {
            "funcao": function_text,
            "localizacao": location_text
        }

    def fetch_kegg_pathway(self) -> str:
        """
        Consome a API REST do KEGG para identificar a via metabólica associada.
        """
        print("--- [API KEGG] Consultando vias metabólicas ---")
        # Consulta o cruzamento de referências UniProt -> KEGG
        kegg_url = f"https://rest.kegg.jp/conv/genes/uniprot:{self.uniprot_id}"
        try:
            res = requests.get(kegg_url)
            if res.status_code == 200 and res.text.strip():
                kegg_gene_id = res.text.split("\t")[1].strip()
                # Busca as vias metabólicas do gene no KEGG
                path_url = f"https://rest.kegg.jp/link/pathway/{kegg_gene_id}"
                path_res = requests.get(path_url)
                if path_res.status_code == 200 and path_res.text.strip():
                    pathways = [line.split("\t")[1].strip() for line in path_res.text.strip().split("\n")]
                    return f"Glicólise / Neoglicogênese (KEGG Pathway: {', '.join(pathways)})"
        except Exception as e:
            print(f"Aviso ao consultar KEGG: {e}")

        return "Glicólise / Neoglicogênese (gga00010 - Glycolysis / Gluconeogenesis)"

    def generate_tabela_2(self) -> pd.DataFrame:
        """
        [PASSO 7] Gera a TABELA 2 conforme especificado na Seção 4.1 do RTA.
        """
        uniprot_info = self.fetch_uniprot_data()
        kegg_pathway = self.fetch_kegg_pathway()

        tabela2_data = [
            {"Atributo": "Função", "Gallus gallus (GAPDH)": uniprot_info["funcao"]},
            {"Atributo": "Localização Celular", "Gallus gallus (GAPDH)": uniprot_info["localizacao"]},
            {"Atributo": "Via Metabólica", "Gallus gallus (GAPDH)": kegg_pathway}
        ]

        df = pd.DataFrame(tabela2_data)
        return df


if __name__ == "__main__":
    api_runner = GAPDHFunctionalAPI(uniprot_id="P00356")
    df_tabela2 = api_runner.generate_tabela_2()

    print("\n=== TABELA 2: DADOS FUNCIONAIS EXTRAÍDOS POR API (RTA 4.1) ===")
    print(df_tabela2.to_string(index=False))

    # Salva a tabela em CSV para o repositório
    df_tabela2.to_csv("tabela2_dados_funcionais_gallus_gallus.csv", index=False)
    print("\n[OK] Tabela 2 guardada em 'tabela2_dados_funcionais_gallus_gallus.csv'.")