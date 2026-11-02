# Trigonometric Least Squares (TLS)

Tiểu luận môn Giải tích số: xấp xỉ dữ liệu có chu kỳ bằng tổng các hàm sin và cos, giải bằng bình phương tối thiểu.

$$\hat y(x) = a_0 + \sum_{k=1}^{K} \left[ a_k \cos(k\omega x) + b_k \sin(k\omega x) \right], \qquad \omega = \frac{2\pi}{T}$$

## Nội dung

| File | Mô tả |
|---|---|
| `tls_synthetic.ipynb` | Kiểm tra TLS trên dữ liệu giả lập, nơi biết trước đáp án đúng |
| `tls_experiment.ipynb` | Áp dụng TLS cho nhiệt độ trung bình ngày tại Quy Nhơn, 01/2016 – 10/2026 |
| `fetch_data.py` | Tải dữ liệu ERA5 qua Open-Meteo vào thư mục `data/` |

## Chạy lại

```bash
pip install numpy pandas matplotlib seaborn jupyter
python fetch_data.py   # không bắt buộc, dữ liệu đã có sẵn trong data/
```

## Kết quả chính

Khớp trên 2016 – 2024, kiểm tra trên 01/2025 – 10/2026, với $T = 365.25$ ngày và $K = 4$:

| Mô hình | RMSE kiểm tra (°C) |
|---|---|
| Hằng số (trung bình tập huấn luyện) | 2.666 |
| Trung bình khí hậu theo ngày | 1.102 |
| **TLS** | **1.022** |

## Nguồn dữ liệu

ERA5 ([Hersbach et al., 2020](https://doi.org/10.1002/qj.3803)), lấy qua [Open-Meteo](https://open-meteo.com/) (CC BY 4.0).
