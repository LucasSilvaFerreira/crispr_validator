#!/usr/bin/env python3

from __future__ import annotations

import base64
import csv
import hashlib
import json
import os
import shutil
import ssl
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence
from urllib.error import HTTPError, URLError
from urllib.parse import urljoin, urlparse
from urllib.request import HTTPRedirectHandler, HTTPSHandler, Request, build_opener, urlopen

try:
    import certifi
except ImportError:  # pragma: no cover - optional dependency
    certifi = None


PORTAL_BASE = "https://api.data.igvf.org"
ACTIVE_FILE_STATUSES = {"in progress", "preview", "released"}
INACTIVE_FILE_STATUSES = {"archived", "deleted", "replaced", "revoked"}

FILE_SET_TYPE_TO_MODALITY = {
    "experimental data": "scRNA",
    "grna sequencing": "gRNA",
    "cell hashing barcode sequencing": "hash",
}

FALLBACK_MODALITY_MAP = {
    "scrna": "rna_seqspec",
    "grna": "sgrna_seqspec",
    "hash": "hash_seqspec",
}

SSL_CONTEXT = ssl.create_default_context(cafile=certifi.where()) if certifi else None
REDIRECT_HTTP_CODES = {301, 302, 303, 307, 308}


@dataclass(frozen=True)
class IGVFCredentials:
    key: str
    secret: str


@dataclass(frozen=True)
class PortalFileReference:
    accession: str
    url: str
    md5sum: str
    filename: str


class _NoRedirectHandler(HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        return None


def build_auth_headers(credentials: Optional[IGVFCredentials]) -> Dict[str, str]:
    headers = {
        "Accept": "application/json",
        "User-Agent": "crispr-validator/igvf",
    }
    if credentials:
        token = base64.b64encode(f"{credentials.key}:{credentials.secret}".encode("utf-8")).decode("ascii")
        headers["Authorization"] = f"Basic {token}"
    return headers


def load_credentials(
    keypair_path: Optional[str] = None,
    api_key: Optional[str] = None,
    api_secret: Optional[str] = None,
    required: bool = False,
) -> Optional[IGVFCredentials]:
    if api_key or api_secret:
        if not (api_key and api_secret):
            raise ValueError("Both IGVF API key and secret are required when specifying credentials explicitly.")
        return IGVFCredentials(api_key, api_secret)
    if keypair_path:
        with open(keypair_path, "r", encoding="utf-8") as handle:
            data = json.load(handle)
        key = str(data.get("key", "")).strip()
        secret = str(data.get("secret", "")).strip()
        if not (key and secret):
            raise ValueError(f"Keypair file {keypair_path} must contain 'key' and 'secret'.")
        return IGVFCredentials(key, secret)
    env_key = os.getenv("IGVF_API_KEY", "").strip()
    env_secret = os.getenv("IGVF_SECRET_KEY", "").strip()
    if env_key and env_secret:
        return IGVFCredentials(env_key, env_secret)
    if required:
        raise RuntimeError(
            "No IGVF credentials available. Provide --igvf-api-key/--igvf-api-secret, --igvf-keypair, "
            "or set IGVF_API_KEY and IGVF_SECRET_KEY."
        )
    return None


def portal_json(path_or_url: str, credentials: Optional[IGVFCredentials]) -> Dict[str, object]:
    url = path_or_url if path_or_url.startswith("http") else urljoin(PORTAL_BASE, path_or_url)
    request = Request(url, headers=build_auth_headers(credentials))
    try:
        with urlopen(request, context=SSL_CONTEXT) as response:
            return json.loads(response.read().decode("utf-8"))
    except HTTPError as exc:
        message = exc.read().decode("utf-8", errors="replace")
        raise RuntimeError(f"IGVF API request failed for {url}: {exc.code} {message}") from exc
    except URLError as exc:
        raise RuntimeError(f"Could not reach IGVF API for {url}: {exc.reason}") from exc


def write_samplesheet(rows: Sequence[Dict[str, str]], output_path: str | Path) -> Path:
    if not rows:
        raise ValueError("No rows available to write.")
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return path


def calculate_md5(path: str | Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.md5()
    with open(path, "rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def _open_without_redirect(request: Request):
    handlers = [_NoRedirectHandler()]
    if SSL_CONTEXT:
        handlers.insert(0, HTTPSHandler(context=SSL_CONTEXT))
    return build_opener(*handlers).open(request)


def _download_portal_redirect(
    url: str,
    destination: Path,
    base_headers: Dict[str, str],
    credentials: IGVFCredentials,
) -> None:
    auth_headers = dict(base_headers)
    auth_headers["Authorization"] = build_auth_headers(credentials)["Authorization"]
    request = Request(url, headers=auth_headers)
    try:
        with _open_without_redirect(request) as response, destination.open("wb") as handle:
            shutil.copyfileobj(response, handle)
        return
    except HTTPError as exc:
        if exc.code not in REDIRECT_HTTP_CODES:
            message = exc.read().decode("utf-8", errors="replace")
            raise RuntimeError(f"Failed to download {url}: {exc.code} {message}") from exc
        redirect_target = exc.headers.get("Location", "").strip()
        if not redirect_target:
            raise RuntimeError(f"Failed to download {url}: missing redirect target from portal response.") from exc

    redirected_url = urljoin(url, redirect_target)
    redirected_request = Request(redirected_url, headers=base_headers)
    try:
        with urlopen(redirected_request, context=SSL_CONTEXT) as response, destination.open("wb") as handle:
            shutil.copyfileobj(response, handle)
    except HTTPError as exc:
        message = exc.read().decode("utf-8", errors="replace")
        raise RuntimeError(f"Failed to download {redirected_url}: {exc.code} {message}") from exc
    except URLError as exc:
        raise RuntimeError(f"Could not download {redirected_url}: {exc.reason}") from exc


def download_file(
    url: str,
    destination: Path,
    credentials: Optional[IGVFCredentials],
    *,
    byte_range_end: Optional[int] = None,
    expected_md5: Optional[str] = None,
) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    base_headers = {
        "User-Agent": "crispr-validator/igvf",
        "Accept": "*/*",
    }
    if byte_range_end is not None and byte_range_end >= 0:
        base_headers["Range"] = f"bytes=0-{byte_range_end}"

    if credentials and url.startswith(PORTAL_BASE):
        _download_portal_redirect(url, destination, base_headers, credentials)
        if expected_md5 and byte_range_end is None:
            actual = calculate_md5(destination)
            if actual != expected_md5:
                raise RuntimeError(
                    f"Checksum mismatch for {destination}. Expected {expected_md5}, observed {actual}."
                )
        return

    last_error: Optional[Exception] = None
    request_headers = [base_headers]
    if credentials:
        request_headers.append({**base_headers, **build_auth_headers(credentials)})

    for headers in request_headers:
        request = Request(url, headers=headers)
        try:
            with urlopen(request, context=SSL_CONTEXT) as response, destination.open("wb") as handle:
                shutil.copyfileobj(response, handle)
            break
        except HTTPError as exc:
            message = exc.read().decode("utf-8", errors="replace")
            last_error = RuntimeError(f"Failed to download {url}: {exc.code} {message}")
            if exc.code not in {401, 403}:
                raise last_error from exc
            if headers is request_headers[-1]:
                raise last_error from exc
        except URLError as exc:
            last_error = RuntimeError(f"Could not download {url}: {exc.reason}")
            if headers is request_headers[-1]:
                raise last_error from exc
    else:  # pragma: no cover - defensive guard
        if last_error:
            raise last_error
        raise RuntimeError(f"Failed to download {url}")

    if expected_md5 and byte_range_end is None:
        actual = calculate_md5(destination)
        if actual != expected_md5:
            raise RuntimeError(
                f"Checksum mismatch for {destination}. Expected {expected_md5}, observed {actual}."
            )


def canonical_modality(file_set_type: str) -> Optional[str]:
    return FILE_SET_TYPE_TO_MODALITY.get(file_set_type.strip().lower())


def extract_accession(value: str) -> str:
    return value.strip("/").split("/")[-1]


def first_matching_file(
    links: Iterable[str],
    credentials: IGVFCredentials,
    *,
    required: bool = True,
) -> Optional[Dict[str, object]]:
    for link in links:
        obj = portal_json(f"{link}@@object?format=json", credentials)
        status = str(obj.get("status", "")).strip().lower()
        if status in INACTIVE_FILE_STATUSES:
            continue
        return obj
    if required:
        raise ValueError("No active IGVF file object matched the requested links.")
    return None


def build_download_url(file_object: Dict[str, object]) -> str:
    href = str(file_object.get("href", "")).strip()
    accession = str(file_object.get("accession", "")).strip()
    object_path = str(file_object.get("@id", "")).strip()
    filename = str(file_object.get("submitted_file_name", "")).strip()
    basename = Path(filename).name if filename else accession
    if href:
        href_name = Path(urlparse(href).path).name
        if "." in href_name or not basename:
            return urljoin(PORTAL_BASE, href)
        if not object_path:
            raise ValueError(f"Could not derive download URL for IGVF file object: {file_object}")
        return urljoin(PORTAL_BASE, f"{object_path}@@download/{basename}")
    if not accession or not object_path:
        raise ValueError(f"Could not derive download URL for IGVF file object: {file_object}")
    return urljoin(PORTAL_BASE, f"{object_path}@@download/{basename}")


def build_file_reference(file_object: Dict[str, object]) -> PortalFileReference:
    accession = str(file_object.get("accession", "")).strip()
    if not accession:
        raise ValueError(f"Missing accession on IGVF file object: {file_object}")
    url = build_download_url(file_object)
    filename = Path(urlparse(url).path).name
    md5sum = str(file_object.get("content_md5sum") or file_object.get("md5sum") or "").strip()
    return PortalFileReference(
        accession=accession,
        url=url,
        md5sum=md5sum,
        filename=filename,
    )


def fallback_seqspec_path(
    modality: str,
    *,
    hash_seqspec: Optional[str],
    rna_seqspec: Optional[str],
    sgrna_seqspec: Optional[str],
) -> str:
    key = FALLBACK_MODALITY_MAP.get(modality.strip().lower())
    options = {
        "hash_seqspec": hash_seqspec or "",
        "rna_seqspec": rna_seqspec or "",
        "sgrna_seqspec": sgrna_seqspec or "",
    }
    return options.get(key, "")


def select_valid_seqspec(
    links: Sequence[str],
    credentials: IGVFCredentials,
) -> Optional[PortalFileReference]:
    for link in links:
        obj = portal_json(f"{link}@@object?format=json", credentials)
        status = str(obj.get("status", "")).strip().lower()
        upload_status = str(obj.get("upload_status", "")).strip().lower()
        if status not in ACTIVE_FILE_STATUSES or upload_status != "validated":
            continue
        return build_file_reference(obj)
    return None


def generate_samplesheet_rows(
    analysis_set_accession: str,
    credentials: IGVFCredentials,
    *,
    hash_seqspec: Optional[str] = None,
    rna_seqspec: Optional[str] = None,
    sgrna_seqspec: Optional[str] = None,
) -> List[Dict[str, str]]:
    analysis_set = portal_json(f"/analysis-sets/{analysis_set_accession}/@@object?format=json", credentials)
    input_file_sets = [
        str(link)
        for link in analysis_set.get("input_file_sets", [])
        if str(link).startswith("/measurement-sets/") or str(link).startswith("/auxiliary-sets/")
    ]
    if not input_file_sets:
        raise ValueError(f"Analysis set {analysis_set_accession} did not expose measurement/auxiliary input file sets.")

    construct_library_sets = [str(link) for link in analysis_set.get("construct_library_sets", [])]
    if len(construct_library_sets) != 1:
        raise ValueError(
            f"Analysis set {analysis_set_accession} must resolve to exactly one construct library set; "
            f"found {len(construct_library_sets)}."
        )
    construct_library = portal_json(f"{construct_library_sets[0]}@@embedded?format=json", credentials)
    guide_file = None
    for file_object in construct_library.get("integrated_content_files", []):
        status = str(file_object.get("status", "")).strip().lower()
        if (
            str(file_object.get("content_type", "")).strip().lower() == "guide rna sequences"
            and status in ACTIVE_FILE_STATUSES
        ):
            file_link = str(file_object.get("@id", "")).strip()
            guide_source = portal_json(f"{file_link}@@object?format=json", credentials) if file_link else file_object
            guide_file = build_file_reference(guide_source)
            break
    if guide_file is None:
        raise ValueError(f"Analysis set {analysis_set_accession} is missing an active guide RNA sequences file.")

    measurement_links = [link for link in input_file_sets if link.startswith("/measurement-sets/")]
    auxiliary_links = [link for link in input_file_sets if link.startswith("/auxiliary-sets/")]

    measurement_properties: Dict[str, Dict[str, object]] = {}
    per_sample_rows: List[Dict[str, str]] = []

    for file_set_link in [*measurement_links, *auxiliary_links]:
        file_set = portal_json(f"{file_set_link}@@object?format=json", credentials)
        file_set_type = str(file_set.get("file_set_type", "")).strip()
        modality = canonical_modality(file_set_type)
        if not modality:
            continue

        measurement_sets = ""
        barcode_onlist_links = []
        onlist_method = ""
        strand_specificity = ""
        if file_set_link.startswith("/measurement-sets/"):
            measurement_sets = str(file_set.get("accession", "")).strip()
            barcode_onlist_links = [str(link) for link in file_set.get("onlist_files", [])]
            onlist_method = str(file_set.get("onlist_method", "")).strip()
            strand_specificity = str(file_set.get("strand_specificity", "")).strip()
            measurement_properties[file_set_link] = {
                "barcode_onlist_links": barcode_onlist_links,
                "onlist_method": onlist_method,
                "strand_specificity": strand_specificity,
            }
        else:
            linked_measurements = [
                str(link) for link in file_set.get("measurement_sets", []) if str(link) in measurement_properties
            ]
            if not linked_measurements:
                raise ValueError(f"Auxiliary set {file_set_link} did not map to a measurement set in this analysis set.")
            props = measurement_properties[linked_measurements[0]]
            barcode_onlist_links = list(props["barcode_onlist_links"])
            onlist_method = str(props["onlist_method"])
            strand_specificity = str(props["strand_specificity"])
            measurement_sets = ", ".join(extract_accession(link) for link in linked_measurements)

        if not barcode_onlist_links:
            raise ValueError(f"Missing barcode onlist on {file_set_link}.")
        if not strand_specificity:
            raise ValueError(f"Missing strand_specificity on {file_set_link}.")
        if not onlist_method:
            raise ValueError(f"Missing onlist_method on {file_set_link}.")
        if onlist_method != "no combination" and len(barcode_onlist_links) != 1:
            raise ValueError(
                f"Unsupported onlist_method '{onlist_method}' on {file_set_link} with "
                f"{len(barcode_onlist_links)} onlist files. The validator requires exactly one resolved barcode onlist."
            )

        barcode_onlist = build_file_reference(first_matching_file(barcode_onlist_links, credentials))

        hash_map = None
        if modality == "hash":
            barcode_map_link = str(file_set.get("barcode_map", "")).strip()
            if not barcode_map_link:
                raise ValueError(f"Hash modality on {file_set_link} is missing barcode_map.")
            hash_map = build_file_reference(portal_json(f"{barcode_map_link}@@object?format=json", credentials))

        reads_by_key: Dict[tuple[str, str, str, str], Dict[str, Dict[str, object]]] = {}
        for file_link in file_set.get("files", []):
            file_object = portal_json(f"{file_link}@@object?format=json", credentials)
            if not str(file_object.get("@id", "")).startswith("/sequence-files/"):
                continue
            if str(file_object.get("status", "")).strip().lower() in INACTIVE_FILE_STATUSES:
                continue
            if str(file_object.get("content_type", "")).strip().lower() != "reads":
                continue
            read_type = str(file_object.get("illumina_read_type", "")).strip().upper()
            if read_type not in {"R1", "R2"}:
                continue
            key = (
                str(file_object.get("sequencing_run", "")).strip(),
                str(file_object.get("lane", "")).strip(),
                str(file_object.get("flowcell_id", "")).strip(),
                str(file_object.get("index", "")).strip(),
            )
            reads_by_key.setdefault(key, {})[read_type] = file_object

        for (sequencing_run, lane, flowcell_id, index), reads in sorted(reads_by_key.items()):
            read1_object = reads.get("R1")
            read2_object = reads.get("R2")
            if not read1_object or not read2_object:
                continue

            seqspec_ref = select_valid_seqspec(
                [str(link) for link in read1_object.get("seqspecs", [])],
                credentials,
            )
            seqspec_path = ""
            seqspec_accession = ""
            seqspec_md5sum = ""
            if seqspec_ref:
                seqspec_path = seqspec_ref.url
                seqspec_accession = seqspec_ref.accession
                seqspec_md5sum = seqspec_ref.md5sum
            else:
                seqspec_path = fallback_seqspec_path(
                    modality,
                    hash_seqspec=hash_seqspec,
                    rna_seqspec=rna_seqspec,
                    sgrna_seqspec=sgrna_seqspec,
                )
            if not seqspec_path:
                raise ValueError(
                    f"Missing seqspec for modality {modality} in analysis set {analysis_set_accession}; "
                    "provide a fallback seqspec path."
                )

            read1 = build_file_reference(read1_object)
            read2 = build_file_reference(read2_object)
            row = {
                "R1_path": read1.url,
                "R1_accession": read1.accession,
                "R1_md5sum": read1.md5sum,
                "R2_path": read2.url,
                "R2_accession": read2.accession,
                "R2_md5sum": read2.md5sum,
                "file_modality": modality,
                "file_set": str(file_set.get("accession", "")).strip(),
                "file_set_type": file_set_type,
                "measurement_sets": measurement_sets,
                "sequencing_run": sequencing_run,
                "lane": lane,
                "flowcell_id": flowcell_id,
                "index": index,
                "seqspec": seqspec_path,
                "seqspec_accession": seqspec_accession,
                "seqspec_md5sum": seqspec_md5sum,
                "barcode_onlist": barcode_onlist.url,
                "barcode_onlist_accession": barcode_onlist.accession,
                "barcode_onlist_md5sum": barcode_onlist.md5sum,
                "onlist_method": onlist_method,
                "strand_specificity": strand_specificity,
                "guide_design": guide_file.url,
                "guide_design_accession": guide_file.accession,
                "guide_design_md5sum": guide_file.md5sum,
                "barcode_hashtag_map": hash_map.url if hash_map else "",
                "barcode_hashtag_map_accession": hash_map.accession if hash_map else "",
                "barcode_hashtag_map_md5sum": hash_map.md5sum if hash_map else "",
            }
            per_sample_rows.append(row)

    if not per_sample_rows:
        raise ValueError(f"No supported per-sample rows could be derived for analysis set {analysis_set_accession}.")
    return per_sample_rows
