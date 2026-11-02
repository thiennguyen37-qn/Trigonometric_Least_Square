"""Tải dữ liệu nhiệt độ trung bình ngày tại Quy Nhơn (ERA5) từ Open-Meteo.

Chạy:     py fetch_data.py
Kết quả:  data/quy_nhon_era5_daily.csv   (ngày, nhiệt độ °C)
          data/quy_nhon_era5_daily.json  (nguồn, tham số, thời điểm tải)

Chỉ dùng thư viện chuẩn của Python, không cần cài thêm gói.
"""

import csv
import json
import sys
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

API_URL = "https://archive-api.open-meteo.com/v1/archive"
PARAMS = {
    "latitude": 13.80,               # Quy Nhơn
    "longitude": 109.20,
    "start_date": "2016-01-01",
    "end_date": "2026-10-25",        # cố định để tái lập được kết quả
    "daily": "temperature_2m_mean",
    "models": "era5",
    "timezone": "Asia/Bangkok",      # GMT+7
}

DATA_DIR = Path(__file__).resolve().parent / "data"
CSV_PATH = DATA_DIR / "quy_nhon_era5_daily.csv"
META_PATH = DATA_DIR / "quy_nhon_era5_daily.json"


def fetch(params):
    """Gọi API và trả về (url, dữ liệu JSON)."""
    url = f"{API_URL}?{urllib.parse.urlencode(params)}"
    try:
        with urllib.request.urlopen(url, timeout=60) as response:
            payload = json.load(response)
    except urllib.error.HTTPError as err:
        body = err.read().decode(errors="replace")
        raise RuntimeError(f"Open-Meteo trả về lỗi {err.code}: {body}") from err
    return url, payload


def main():
    sys.stdout.reconfigure(encoding="utf-8")  # in được tiếng Việt trên Windows
    url, payload = fetch(PARAMS)
    dates = payload["daily"]["time"]
    temps = payload["daily"]["temperature_2m_mean"]

    DATA_DIR.mkdir(exist_ok=True)
    with CSV_PATH.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["date", "temperature_2m_mean"])
        writer.writerows(zip(dates, temps))

    meta = {
        "source": "Open-Meteo Historical Weather API (ERA5 reanalysis)",
        "url": url,
        "downloaded_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "grid_latitude": payload["latitude"],
        "grid_longitude": payload["longitude"],
        "elevation_m": payload["elevation"],
        "timezone": payload["timezone"],
        "units": payload["daily_units"],
        "n_days": len(dates),
        "n_missing": sum(t is None for t in temps),
    }
    META_PATH.write_text(
        json.dumps(meta, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    print(f"Đã lưu {meta['n_days']} ngày ({dates[0]} → {dates[-1]}) vào {CSV_PATH}")
    print(f"Ô lưới ERA5: {meta['grid_latitude']}°N, {meta['grid_longitude']}°E")
    print(f"Số giá trị thiếu: {meta['n_missing']}")


if __name__ == "__main__":
    main()
