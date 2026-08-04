import csv
import logging
import requests
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional
import concurrent.futures


@dataclass
class SraMetadata:
    ids: List[str] = field(default_factory=list)
    mems: List[str] = field(default_factory=list)
    real_mems: List[str] = field(default_factory=list)

    @classmethod
    def from_tsv(cls, tsv_file: str | Path) -> "SraMetadata":
        """
        Parses a TSV file and initializes the dataclass.
        Skips the header and strips the '.tsv' extension from the ID column.
        """
        metadata = cls()
        with open(tsv_file, mode="r", encoding="utf-8") as f:
            reader = csv.reader(f, delimiter="\t")
            next(reader) 
            
            for row in reader:
                if len(row) >= 3:
                    srr_id = row[0].replace(".tsv", "").strip()
                    metadata.ids.append(srr_id)
                    metadata.mems.append(row[1].strip())
                    metadata.real_mems.append(row[2].strip())
        return metadata

@dataclass
class FastqDownloader:
    """Handles downloading FASTQ files from the ENA API based on SRA metadata."""
    
    metadata : SraMetadata
    output_dir : Path = Path("database/raw_data").mkdir(parents=True, exist_ok=True)
    ena_api_url: str = "https://www.ebi.ac.uk/ena/portal/api/filereport"  
    max_workers: int = 4 

    def fetch_data(self, max_bytes: Optional[int] = None) -> None:
        """
        Iterates through the IDs and downloads their corresponding FASTQ files.
        
        Args:
            max_bytes: If set, only downloads the first N bytes for testing purposes.
        """
        logging.info(f"Starting parallel download with {self.max_workers} workers...")
        with concurrent.futures.ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            futures = {
                executor.submit(self._process_single_id, id, max_bytes): id for id in self.metadata.ids
            }

            for future in concurrent.futures.as_completed(futures):
                id = futures[future]
                try:
                    future.result()
                except Exception as e:
                    logging.error(f"[{id}] A critical error occurred: {e}")


    def _process_single_id(self, id: str, max_bytes: Optional[int] = None) -> None:
        logging.debug(f"[{id}] Fetching download links...")
        urls = self._get_ena_fastq_urls(id)  
        if not urls:
            logging.warning(f"[{id}] No FASTQ URLs found.")
            return

        for url in urls:
            self._download_file(url, max_bytes=max_bytes)

    def _get_ena_fastq_urls(self, id: str) -> List[str]:
        """
        Queries the ENA API to get direct FTP/HTTP download links for FASTQ files.
        """
        params = {
            "accession": id,
            "result": "read_run",
            "fields": "fastq_ftp",
            "format": "json"
        }
        
        try:
            response = requests.get(self.ena_api_url, params=params, timeout=10)
            response.raise_for_status()
            data = response.json()
            
            if not data:
                return []

            ftp_links_str = data[0].get("fastq_ftp", "")
            if not ftp_links_str:
                return []

            return [f"http://{link}" for link in ftp_links_str.split(";")]
            
        except requests.RequestException as e:
            logging.error(f"Error fetching metadata for {id}: {e}")
            return []

    def _download_file(self, url: str, max_bytes: Optional[int] = None) -> None:
        filename = url.split("/")[-1]
        
        if max_bytes:
            filename = f"partial_{filename}"
            
        output_path = self.output_dir / filename

        if output_path.exists():
            logging.debug(f"File {filename} already exists. Skipping download.")
            return

        logging.debug(f"Downloading {filename} from {url}...")
        
        try:
            bytes_downloaded = 0
            with requests.get(url, stream=True, timeout=30) as response:
                response.raise_for_status()
                
                with open(output_path, "wb") as f:
                    for chunk in response.iter_content(chunk_size=8192):
                        if chunk:
                            f.write(chunk)
                            bytes_downloaded += len(chunk)
                            
                            if max_bytes and bytes_downloaded >= max_bytes:
                                logging.debug(f"Reached test limit of {max_bytes} bytes. Stopping download.")
                                break
                                
            logging.info(f"Successfully downloaded {filename} ({bytes_downloaded} bytes)")
            
        except requests.RequestException as e:
            logging.error(f"Failed to download {filename}: {e}")
            if output_path.exists():
                output_path.unlink()


if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG, format='%(levelname)s: %(message)s')
    
    script_dir = Path(__file__).resolve().parent
    project_root = script_dir.parent.parent 
    database_dir = project_root / "database"
    metadata_path = database_dir / "metadata.tsv"
    
    try:
        metadata = SraMetadata.from_tsv(metadata_path)
        
        downloader = FastqDownloader(metadata=metadata, output_dir=database_dir / "raw_data", max_workers=4)
        
        # TEST MODE: Set max_bytes to 1 MB (1024 * 1024 bytes)
        logging.info("Starting download in TEST MODE (1MB limit per file)")
        downloader.fetch_data(max_bytes=1048576) 
        
    except FileNotFoundError:
        logging.error(f"Could not find {metadata_path}. Please create it first.")