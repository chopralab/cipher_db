import json
import os
from pprint import pprint
import requests
import time
from typing import Dict, List, Optional
import urllib3


urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

ASKCOS_HOST = "https://35.203.15.8"
USERNAME = "cipher"
PASSWORD = "password"

# TODO(degraff): these environment variables will need to be set properly on the machine
os.environ["ASKCOS_HOST"] = ASKCOS_HOST
os.environ["ASKCOS_USERNAME"] = USERNAME
os.environ["ASKCOS_PASSWORD"] = PASSWORD


class AskcosClient:
    AUTHENTICATION_ENDPOINT = "/api/v2/token-auth/"
    TREE_BUILDER_ENDPOINT = "/api/v2/tree-builder/"
    SC_SCORE_ENDPOINT = "/api/v2/scscore/"
    TASK_RETRIEVAL_ENDPOINT = "/api/v2/celery/task/"

    def __init__(
        self,
        host: Optional[str] = None,
        params: Optional[Dict] = None,
        timeout: int = 600,
        interval: int = 5,
        authenticate: bool = False,
    ):
        self.host = host or os.environ["ASKCOS_HOST"]
        self.params = params.copy() if params is not None else {}
        self.timeout = timeout
        self.interval = interval
        self.headers = {}
        if authenticate:
            self.authenticate()

    def authenticate(self):
        url = self.host + AskcosClient.AUTHENTICATION_ENDPOINT
        resp = requests.post(
            url,
            data={
                "username": os.environ["ASKCOS_USERNAME"],
                "password": os.environ["ASKCOS_PASSWORD"],
            },
            verify=False,
        )
        token = resp.json()["token"]
        self.headers = {"Authorization": f"Bearer {token}"}

    def get_trees(self, smi: str, params: Optional[Dict] = None) -> Optional[List[Dict]]:
        """get the retrosynthetic trees for the given smiles

        Parameters
        ----------
        smi : str
            the SMILES string
        params : Optional[Dict], default=None
            an optional dictionary containing the parameters of the query. If None, use the default
            parameters of this TreeBuilder

        Returns
        -------
        Optional[Dict]
            the JSON dictionaries corresponding to retrosynthetic trees of the tree builder
            request. None if the tree builder job failed for any reason
        """
        url = self.host + AskcosClient.TREE_BUILDER_ENDPOINT
        if params is None:
            params = dict(smiles=smi, **self.params)
        else:
            params = dict(smiles=smi, **params)

        try:
            resp = requests.post(url, data=params, timeout=20, verify=False)
        except requests.Timeout:
            print(f"Error submitting tree job for SMILES: {smi}")
            return None

        task_id = resp.json()["task_id"]

        return self.get_result(task_id)

    def get_result(self, task_id) -> Optional[List[Dict]]:
        url = self.host + AskcosClient.TASK_RETRIEVAL_ENDPOINT + f"{task_id}/"

        start = time.time()
        while True:
            resp = requests.get(url, verify=False)
            res = resp.json()

            if res["complete"]:
                return res["output"]
            if res["failed"] or ((time.time() - start) > self.timeout):
                return None

            time.sleep(self.interval)

    def sc_score(self, smi: str) -> float:
        url = self.host + AskcosClient.SC_SCORE_ENDPOINT
        try:
            resp = requests.post(url, data={"smiles": smi}, timeout=20, verify=False)
        except requests.Timeout:
            print(f"Error submitting tree job for SMILES: {smi}")
            return None

        return resp.json()["score"]


if __name__ == "__main__":
    params = {"version": 1, "expansion_time": 60}
    client = AskcosClient(None, params)
    smi = "CCN(CC)CCOC(C)(c1ccccc1)c1ccc(Cl)cc1"

    trees = client.get_trees(smi)
    print(json.dumps(trees, indent=2))

    print(client.sc_score(smi))
    exit()
