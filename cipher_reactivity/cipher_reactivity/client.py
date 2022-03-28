from enum import Enum
import json
import os
import time
from typing import Any, Dict, Iterable, List, Optional

import requests
from requests.auth import AuthBase
from requests_toolbelt import sessions
import urllib3

urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)


class AskcosEndpoints(Enum):
    AUTHENTICATION = "/api/v2/token-auth/"
    FORWARD = "/api/v2/forward/"
    SC_SCORE = "/api/v2/scscore/"
    TREE_BUILDER = "/api/v2/tree-builder/"
    TASK_RETRIEVAL = "/api/v2/celery/task/"


class AskcosClient:
    def __init__(
        self,
        host: Optional[str] = None,
        auth: Optional[AuthBase] = None,
        tree_params: Optional[Dict] = None,
        timeout: int = 240,
        interval: int = 5,
        authenticate: bool = False,
    ):
        try:
            self.sess = sessions.BaseUrlSession(host or os.environ["ASKCOS_HOST"])
        except KeyError:
            raise ValueError(
                'arg "host" was not supplied, but ASKCOS_HOST environment variable is not set!'
            )
        self.sess.auth = auth

        if self.sess.get("", verify=False).status_code != 200:
            msg = "ASKCOS host does not exist! "
            if host is not None:
                msg += f"Tried to connect to: {self.sess.base_url}"
            else:
                msg += f"Used value of ASKCOS_HOST environment variable: {self.sess.base_url}"
            raise ValueError(msg)

        self.tree_params = tree_params.copy() if tree_params is not None else {}
        self.timeout = timeout
        self.interval = interval
        
        if authenticate:
            self.authenticate()

    def authenticate(self):
        # url = self.host + AskcosEndpoints.AUTHENTICATION.value
        resp = self.sess.post(
            AskcosEndpoints.AUTHENTICATION.value,
            data={
                "username": os.environ["ASKCOS_USERNAME"],
                "password": os.environ["ASKCOS_PASSWORD"],
            },
            verify=False,
        )
        token = resp.json()["token"]
        self.sess.headers = {"Authorization": f"Bearer {token}"}

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
            request. None if the tree builder job failed for any server error

        Raises
        ------
        ValueError
            if the request failed due to a client error, e.g., invalid SMILES string or bad tree
            parameter
        """
        # url = self.host + AskcosEndpoints.TREE_BUILDER.value
        if tree_params is None:
            payload = dict(smiles=smi, **self.tree_params)
        else:
            payload = dict(smiles=smi, **tree_params)

        try:
            resp = self.sess.post(
                AskcosEndpoints.TREE_BUILDER.value, data=payload, timeout=20, verify=False
            )
            resp.raise_for_status()
        except requests.HTTPError as e:
            if resp.json()["smiles"] != smi:
                raise ValueError(f"Invalid SMILES string supplied! got: {smi}")
            raise ValueError(e)
        except requests.RequestException as e:
            print(f"Error submitting tree job for SMILES: {smi}")
            print(e)
            return None

        return self.get_result(resp.json()["task_id"], self.interval, self.timeout)

    def predict_products(
        self,
        reactants: Iterable[str],
        reagents: Optional[Iterable[str]] = None,
        solvent: Optional[str] = "",
        num_results: int = 100,
        atommap: bool = False,
    ):
        # url = self.host + AskcosEndpoints.FORWARD.value

        payload = {
            "reactants": ".".join(reactants),
            "reagents": ".".join(reagents) if reagents else "",
            "solvent": solvent,
            "num_results": num_results,
            "atommmap": atommap,
        }

        try:
            resp = self.sess.post(AskcosEndpoints.FORWARD.value, payload, timeout=20, verify=False)
            resp.raise_for_status()
        # except requests.HTTPError as e:
        #     if resp.json()["smiles"] != smi:
        #         raise ValueError(f"Invalid SMILES string supplied! got: {smi}")
        #     raise ValueError(e)
        except requests.RequestException as e:
            print(f"Error submitting tree job for SMILES: {smi}")
            print(e)
            return None

        return self.get_result(resp.json()["task_id"], 1, self.timeout)

    def sc_score(self, smi: str) -> float:
        # url = self.host + AskcosEndpoints.SC_SCORE.value

        try:
            resp = self.sess.post(
                AskcosEndpoints.SC_SCORE.value, data={"smiles": smi}, timeout=20, verify=False
            )
            resp.raise_for_status()
        except requests.HTTPError as e:
            if resp.json()["smiles"] != smi:
                raise ValueError(f"Invalid SMILES string supplied! got: {smi}")
            raise ValueError(e)
        except requests.RequestException as e:
            print(f"Error submitting scscore job for SMILES: {smi}")
            print(e)
            return None

        return resp.json()["score"]

    def get_result(self, task_id, interval: float, timeout: float) -> Optional[Any]:
        # url = self.host + AskcosEndpoints.TASK_RETRIEVAL.value + f"{task_id}/"

        start = time.time()
        while True:
            resp = self.sess.get(AskcosEndpoints.TASK_RETRIEVAL.value + f"{task_id}/", verify=False)
            res = resp.json()

            if res["complete"]:
                return res["output"]
            if res["failed"]:
                print(res["error"])
                return None
            if (time.time() - start) > timeout:
                return None

            time.sleep(interval)


if __name__ == "__main__":
    params = {"version": 1, "expansion_time": 60}
    client = AskcosClient(None, params)
    smi = "CCN(CC)CCOC(C)(c1ccccc1)c1ccc(Cl)cc1"

    trees = client.get_trees(smi)
    print(json.dumps(trees, indent=2))

    print(client.sc_score(smi))
    exit()
