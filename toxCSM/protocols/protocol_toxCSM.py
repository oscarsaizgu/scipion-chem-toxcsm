################ VERÓNICA GAMO PAREJO ##############################

# General imports 
import os
import time
import json
import requests

# Specific imports
from toxCSM import Plugin
from pyworkflow.protocol.params import PointerParam, EnumParam, StringParam
from pwem.protocols import EMProtocol
from pwchem.utils import *
from pwchem import Plugin
from pwchem.constants import RDKIT_DIC
from pyworkflow.utils.path import moveFile


class ProtChemToxCSM(EMProtocol):

    """Toxicity prediction of small ligands with toxCSM
    
    AI Generated:

        ProtChemToxCSM - User Manual

        Overview
        --------
        The ProtChemToxCSM protocol performs structure-based toxicity prediction
        for small molecules using the toxCSM web service. It applies graph-based
        signatures and machine learning models to estimate multiple toxicological
        endpoints directly from the chemical structure of a compound.

        This protocol is designed for early-stage drug discovery workflows, where
        rapid toxicity screening is essential to prioritize safe and viable
        candidates before experimental validation.

        Input Requirements
        ------------------
        1. **Molecule Set**:
           - A SetOfSmallMolecules object containing the compounds to analyze.
           - Molecules must include valid structural information (e.g., SDF, MOL).

        2. **Ligand Selection**:
           - A specific ligand name must be provided.
           - The protocol extracts the corresponding molecule from the input set.

        3. **Prediction Type**:
           - Defines the toxicity category to evaluate:
             - stress_response
             - nuclear_response
             - environmental
             - organic
             - human_dose_response
             - genomic
             - all (default; runs all available models)

        Workflow
        --------
        1. **Ligand selection**:
           - The protocol searches the input molecule set for the selected ligand.
           - The corresponding structure file is retrieved.

        2. **SMILES generation**:
           - The molecule is converted into a SMILES representation using RDKit.
           - This representation is required for toxCSM predictions.

        3. **Job submission**:
           - The SMILES string and selected prediction type are sent to the toxCSM API.
           - A unique job ID is returned and stored.

        4. **Job tracking**:
           - The job ID is saved locally for reproducibility and tracking.

        5. **Result retrieval**:
           - After a waiting period, results are retrieved from the toxCSM server.
           - Predictions are returned in JSON format.

        6. **Result storage**:
           - The retrieved results are saved as a JSON file for downstream analysis.

        Outputs
        -------
        - **toxCSM_results.json**:
          - Contains toxicity predictions for the selected ligand.
          - Includes probabilities, classifications, or quantitative values
            depending on the endpoint.

        - **job_id.txt**:
          - Stores the unique identifier of the submitted toxCSM job.
          - Useful for reproducibility and external tracking.

        Validation & Warnings
        ---------------------
        - Ensure that the selected ligand exists in the input molecule set.
        - Molecule files must be structurally valid for SMILES conversion.
        - Internet connection is required to access the toxCSM API.
        - If the API request fails, the protocol raises an exception.
        - The waiting time for results is fixed and may not reflect actual
          server processing time.

        Practical Recommendations
        -------------------------
        - Use the "all" prediction type for a comprehensive toxicity profile.
        - Use specific prediction types to reduce runtime and focus analysis.
        - Verify SMILES correctness if unexpected results occur.
        - Combine toxicity results with docking, ADMET, or activity predictions
          for better compound prioritization.

        Final Perspective
        -----------------
        ProtChemToxCSM provides a fast and automated way to assess the
        toxicological risk of small molecules using structure-based predictions.
        By integrating toxCSM into Scipion-Chem workflows, users can efficiently
        filter unsafe compounds and support decision-making in computational
        drug discovery pipelines.
    """
    
    _label = 'toxicity prediction of small ligands'
    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        self.job_id= None
        self.smi=None


    def _defineParams(self, form):

        form.addSection(label='Input')
                      
        form.addParam('inputSet', PointerParam, pointerClass='SetOfSmallMolecules',
                      label='Molecule for Toxicity Prediction', allowsNull=False,
                      help='Select the set of small molecules for toxicity prediction.')
        
        form.addParam('inputLigand', StringParam, label='Selected just one ligand:',
                      help='Ligand name for toxicity prediction.')
        
        form.addParam('predType', EnumParam, label='Prediction Type', default=6,
                      choices=['stress_response', 'nuclear_response', 'environmental', 'organic', 'human_dose_response', 'genomic', 'all'],
                      help='Type of prediction to perform.')
        
    def _insertAllSteps(self):
        self._insertFunctionStep('extractSmile')
        self._insertFunctionStep('submitJob')
        self._insertFunctionStep('saveJobIdToFile')
        self._insertFunctionStep('retrieveResults')

    
    def extractSmile(self):
        myLigand=str(self.inputLigand)
        for mol in self.inputSet.get():
            if str(mol.getMolName()) in myLigand:
                ligand_file=mol.getFileName()
        if not ligand_file:
            raise ValueError("No ligand provided.")
        
        self.smi = self.getSMI(ligand_file,1)
        print("SMILES for toxicity prediction:", self.smi)
            
    def getSMI(self, mol, nt):
        ''' Generates a SMILES representation of a molecule from a given input file'''
        fnSmall = os.path.abspath(mol)
        fnRoot, ext = os.path.splitext(os.path.basename(fnSmall))

        if ext != '.smi':
            outDir = os.path.abspath(self._getExtraPath())
            fnOut = os.path.abspath(self._getExtraPath(fnRoot + '.smi'))
            args = ' -i "{}" -of smi -o {} --outputDir {} -nt {}'.format(fnSmall, fnOut, outDir, nt)
            Plugin.runScript(self, 'rdkit_IO.py', args, env=RDKIT_DIC, cwd=outDir)    
        return self.parseSMI(fnOut)
        
    def parseSMI(self, smiFile):
        smi = None
        with open(smiFile) as f:
            for line in f:
                smi = line.split()[0].strip()
                if not smi.lower() == 'smiles':
                    break
        return smi
    
    def submitJob(self):
        pred_type_choices = [
            'stress_response', 'nuclear_response', 'environmental',
            'organic', 'human_dose_response', 'genomic', 'all'
        ]

        selected_pred_type = pred_type_choices[self.predType.get()]
        url = "https://biosig.lab.uq.edu.au/toxcsm/api/predict"
        data = {
            'smiles_string': self.smi,
            'pred_type': selected_pred_type
        }

        print(f'Submitting job to {url} with data: {data}')

        try:

            response = requests.post(url, data=data)
            print(f'Response status code: {response.status_code}')
            print(f'Response content: {response.content}')
            response.raise_for_status() 

            job_id = response.json().get('job_id')
            if job_id:
                print('Job submitted successfully. Job ID: {}'.format(job_id))
                self.job_id = job_id
            else:
                print('Job ID not found in the response.')
                print(f'Response content: {response.content}')
                raise Exception('Failed to submit job: Job ID not found')

        except requests.exceptions.RequestException as e:
            print(f'HTTP Request failed: {str(e)}')
            print(f'Response content: {response.content if response else "No response"}')
            raise Exception('Failed to submit job')

        except json.JSONDecodeError as e:
            print(f'Failed to parse JSON response: {str(e)}')
            print(f'Response content: {response.content}')
            raise Exception('Failed to submit job')
        
    def saveJobIdToFile(self):
        filename='job_id.txt'
        if self.job_id is not None:
            with open(filename, 'w') as file:
                file.write(self.job_id)
            fnTxt = self._getExtraPath(filename)
            moveFile(filename, fnTxt)
            print(f"Job ID saved to {filename}")
        else:
            print("No job ID to save.")
        
    def retrieveResults(self):
        result_url = "https://biosig.lab.uq.edu.au/toxcsm/api/predict"
        data = {'job_id': self.job_id}
        time.sleep(120) #Maximum waiting time of te webpage
        response = requests.get(result_url, params=data)
        if response.status_code == 200:
            results = response.json()
            self.saveResultsToJson(results)

    
    def saveResultsToJson(self, parsed_results, filename='toxCSM_results.json'):
        with open(filename, 'w') as json_file:
            json.dump(parsed_results, json_file, indent=4)

        fnJson = self._getExtraPath(filename)
        moveFile(filename, fnJson)

        print(f"Results saved to {filename}")


        



        

    
        

        
        

    
    
