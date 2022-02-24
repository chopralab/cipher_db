from enum import Enum
import json
import os
import time
from typing import Dict, List, Optional
import urllib3

import requests

urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

class AskcosEndpoints(Enum):
    AUTHENTICATION = "/api/v2/token-auth/"
    TREE_BUILDER = "/api/v2/tree-builder/"
    SC_SCORE = "/api/v2/scscore/"
    TASK_RETRIEVAL = "/api/v2/celery/task/"

class AskcosClient:
    def __init__(
        self,
        host: Optional[str] = None,
        tree_params: Optional[Dict] = None,
        timeout: int = 600,
        interval: int = 5,
        authenticate: bool = False,
    ):
        self.host = host or os.environ["ASKCOS_HOST"]
        self.tree_params = tree_params.copy() if tree_params is not None else {}
        self.timeout = timeout
        self.interval = interval
        self.headers = {}
        if authenticate:
            self.authenticate()

    def authenticate(self):
        url = self.host + AskcosEndpoints.AUTHENTICATION.value
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

    def get_trees(self, smi: str, tree_params: Optional[Dict] = None) -> Optional[List[Dict]]:
        """get the retrosynthetic trees for the given smiles

        Parameters
        ----------
        smi : str
            the SMILES string
        tree_params : Optional[Dict], default=None
            an optional dictionary containing the parameters of the query. If None, use the default
            parameters of this TreeBuilder

        Returns
        -------
        Optional[Dict]
            the dictionaries corresponding to retrosynthetic trees of the tree builder
            request. None if the tree builder job failed for any reason
        """
        url = self.host + AskcosEndpoints.TREE_BUILDER.value
        if tree_params is None:
            tree_params = dict(smiles=smi, **self.tree_params)
        else:
            tree_params = dict(smiles=smi, **tree_params)

        try:
            resp = requests.post(url, data=tree_params, timeout=20, verify=False)
        except requests.Timeout:
            print(f"Error submitting tree job for SMILES: {smi}")
            return None

        task_id = resp.json()["task_id"]

        return self.get_result(task_id)

    def get_result(self, task_id) -> Optional[List[Dict]]:
        url = self.host + AskcosEndpoints.TASK_RETRIEVAL.value + f"{task_id}/"

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
        url = self.host + AskcosEndpoints.SC_SCORE.value
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
