# Phương pháp bình phương tối thiểu lượng giác (Trigonometric Least Squares – TLS)

**Tiểu luận cuối kỳ — Giải tích số (Numerical Analysis)**

> **Hướng dẫn sử dụng:** Đây là khung sườn để phát triển thành tiểu luận hoàn chỉnh.
> - Các trường dạng `[...]` là thông tin cần điền trực tiếp.
> - Các khối trích dẫn bắt đầu bằng **Cần bổ sung** là gợi ý nội dung; hãy thay bằng nội dung, công thức, hình ảnh, kết quả thực nghiệm và tài liệu tham khảo của bạn, rồi xóa phần gợi ý.
> - Có thể điều chỉnh tên mục theo yêu cầu của giảng viên.
>
> **Lưu ý ký hiệu:** trong tiểu luận này, TLS là *Trigonometric Least Squares*, không phải *Total Least Squares*.

---

## Thông tin chung

| Mục | Nội dung |
|---|---|
| Tên đề tài | Phương pháp bình phương tối thiểu lượng giác (Trigonometric Least Squares) |
| Học phần | Giải tích số (Numerical Analysis) |
| Sinh viên | [Họ và tên] |
| Mã sinh viên | [Mã sinh viên] |
| Giảng viên | [Tên giảng viên] |
| Ngày hoàn thành | [Ngày/tháng/năm] |

## Tóm tắt

> **Cần bổ sung:** Trình bày ngắn gọn bài toán xấp xỉ dữ liệu, ý tưởng của TLS, phương pháp tính toán được triển khai, thiết kế thực nghiệm Python và những kết quả chính. Viết phần này sau khi hoàn thành các phần còn lại.

**Từ khóa:** Trigonometric Least Squares, bình phương tối thiểu, Fourier, xấp xỉ hàm, sai số số học, Python.

## Mục lục

- **Phần I. Cơ sở lý thuyết**
  1. Giới thiệu bài toán
  2. Ôn tập phương pháp bình phương tối thiểu
  3. Cơ sở toán học của Trigonometric Least Squares
  4. Thuật toán và độ phức tạp tính toán
- **Phần II. Thực nghiệm bằng Python**
  5. Mục tiêu và thiết kế thực nghiệm
  6. Dữ liệu thực nghiệm
  7. Cài đặt TLS bằng Python
  8. Chỉ số đánh giá
  9. Kết quả và thảo luận
  10. Kết luận
- Tài liệu tham khảo
- Phụ lục A. Mã nguồn
- Phụ lục B. Bảng kết quả bổ sung

> Nếu xuất sang Word/PDF, có thể thay mục lục này bằng mục lục tự động của trình soạn thảo.

---

## PHẦN I. CƠ SỞ LÝ THUYẾT

### 1. Giới thiệu bài toán

#### 1.1. Bối cảnh và động cơ

- Giới thiệu bài toán xấp xỉ một tập dữ liệu rời rạc bằng một hàm toán học.
- Nêu những trường hợp dữ liệu có tính chu kỳ hoặc dao động mà mô hình đa thức đơn giản có thể không mô tả thuận tiện.
- Giới thiệu ý tưởng dùng các hàm sin và cos làm hàm cơ sở.

#### 1.2. Mục tiêu của tiểu luận

- Trình bày cơ sở toán học của phương pháp Trigonometric Least Squares.
- Thiết lập bài toán dưới dạng tối ưu hóa tổng bình phương sai số.
- Xây dựng thuật toán tính toán bằng Python.
- Đánh giá kết quả thông qua các bộ dữ liệu và chỉ số sai số phù hợp.

#### 1.3. Phạm vi và cấu trúc

> **Cần bổ sung:** Nêu rõ tiểu luận tập trung vào TLS với chu kỳ đã biết hay có bước ước lượng chu kỳ; số điều hòa được khảo sát; loại dữ liệu và giới hạn của thực nghiệm.

### 2. Ôn tập phương pháp bình phương tối thiểu

#### 2.1. Bài toán xấp xỉ dữ liệu

Cho tập dữ liệu:

$$
(x_i, y_i), \quad i = 1, \ldots, n.
$$

Mục tiêu là xây dựng hàm xấp xỉ $\hat y(x)$ sao cho giá trị dự đoán gần với dữ liệu quan sát.

#### 2.2. Nguyên lý bình phương tối thiểu

Phần dư tại quan sát thứ $i$ được định nghĩa là:

$$
r_i = y_i - \hat y(x_i).
$$

Hàm mục tiêu là tổng bình phương phần dư:

$$
S = \sum_{i=1}^{n} \left[ y_i - \hat y(x_i) \right]^2.
$$

Bài toán bình phương tối thiểu tìm các tham số làm $S$ nhỏ nhất.

#### 2.3. Dạng ma trận của bài toán

Với mô hình tuyến tính theo tham số:

$$
\mathbf{y} \approx A \boldsymbol{\beta},
$$

bài toán có dạng:

$$
\min_{\boldsymbol{\beta}} \left\| \mathbf{y} - A \boldsymbol{\beta} \right\|_2^2.
$$

Nếu ma trận $A$ có hạng cột đầy đủ, nghiệm bình phương tối thiểu là duy nhất và thỏa mãn phương trình chuẩn:

$$
A^{\mathsf T} A \boldsymbol{\beta} = A^{\mathsf T} \mathbf{y}.
$$

> **Cần bổ sung:** Giải thích vì sao trong tính toán thực tế thường ưu tiên các phương pháp QR hoặc SVD hơn việc tính nghịch đảo $(A^{\mathsf T}A)^{-1}$ trực tiếp (số điều kiện của $A^{\mathsf T}A$ bằng bình phương số điều kiện của $A$).

### 3. Cơ sở toán học của Trigonometric Least Squares

#### 3.1. Xây dựng mô hình lượng giác

Với chu kỳ $T$ đã biết, đặt tần số góc:

$$
\omega = \frac{2\pi}{T}.
$$

Mô hình gồm $K$ điều hòa:

$$
\hat y(x) = a_0 + \sum_{k=1}^{K} \left[ a_k \cos(k\omega x) + b_k \sin(k\omega x) \right].
$$

Trong đó:

- $a_0$: thành phần hằng số hoặc mức nền.
- $a_k, b_k$: các hệ số cần ước lượng.
- $T$: chu kỳ cơ bản giả định của dữ liệu.
- $K$: số điều hòa được sử dụng.

> **Cần bổ sung:** Giải thích ý nghĩa của các thành phần và ảnh hưởng của việc tăng $K$.

#### 3.2. Thiết lập hàm mục tiêu

Thay mô hình vào nguyên lý bình phương tối thiểu:

$$
S(a_0, a_1, b_1, \ldots, a_K, b_K) = \sum_{i=1}^{n} \left[ y_i - \hat y(x_i) \right]^2.
$$

Mục tiêu là tìm các hệ số làm $S$ nhỏ nhất.

#### 3.3. Biểu diễn dưới dạng ma trận

Đặt vector tham số:

$$
\boldsymbol{\beta} = [a_0, a_1, b_1, \ldots, a_K, b_K]^{\mathsf T}.
$$

Với mỗi quan sát $x_i$, hàng thứ $i$ của ma trận thiết kế có dạng:

$$
A_i = \left[ 1, \cos(\omega x_i), \sin(\omega x_i), \ldots, \cos(K\omega x_i), \sin(K\omega x_i) \right].
$$

Khi đó:

$$
\mathbf{y} \approx A \boldsymbol{\beta}, \qquad
\hat{\boldsymbol{\beta}} = \arg\min_{\boldsymbol{\beta}} \left\| \mathbf{y} - A \boldsymbol{\beta} \right\|_2^2.
$$

Nếu $A$ có hạng cột đầy đủ, nghiệm bình phương tối thiểu được mô tả bằng phương trình chuẩn:

$$
A^{\mathsf T} A \hat{\boldsymbol{\beta}} = A^{\mathsf T} \mathbf{y}.
$$

> **Cần bổ sung:** Giải thích rằng công thức nghiệm dạng nghịch đảo $\hat{\boldsymbol{\beta}} = (A^{\mathsf T}A)^{-1}A^{\mathsf T}\mathbf{y}$ chỉ áp dụng được khi $A^{\mathsf T}A$ khả nghịch (tức $A$ có hạng cột đầy đủ). Trong Python nên dùng bộ giải bình phương tối thiểu ổn định về số.

#### 3.4. Điều kiện và các vấn đề cần lưu ý

- **Chu kỳ $T$:** Nếu chọn sai chu kỳ, mô hình có thể không phản ánh được cấu trúc thực của dữ liệu.
- **Số điều hòa $K$:** Giá trị nhỏ có thể gây thiếu khớp; giá trị quá lớn có thể khiến mô hình khớp cả nhiễu.
- **Số quan sát và hạng ma trận:** Với $K$ điều hòa có $2K+1$ hệ số, nên cần $n \ge 2K+1$. Tuy nhiên điều kiện này chưa đủ để bảo đảm nghiệm duy nhất nếu các cột của $A$ phụ thuộc tuyến tính.
- **Phân bố điểm lấy mẫu:** Các điểm lấy mẫu và khoảng quan sát ảnh hưởng đến khả năng nhận diện chu kỳ và hệ số.
- **Nhiễu và ngoại lệ:** TLS không tự động loại bỏ nhiễu hoặc ngoại lệ.
- **Ổn định số:** Nên dùng QR/SVD hoặc `numpy.linalg.lstsq` thay vì tự tính `(A.T @ A)^(-1) @ A.T @ y`.

### 4. Thuật toán và độ phức tạp tính toán

#### 4.1. Các bước của thuật toán

1. Nhận dữ liệu $(x_i, y_i)$, chu kỳ $T$ và số điều hòa $K$.
2. Tính $\omega = 2\pi / T$.
3. Xây dựng ma trận thiết kế $A$ từ các cột hằng số, cos và sin.
4. Giải bài toán bình phương tối thiểu để tìm $\hat{\boldsymbol{\beta}}$.
5. Tính giá trị dự đoán $\hat{\mathbf{y}} = A \hat{\boldsymbol{\beta}}$.
6. Tính các chỉ số sai số và trực quan hóa kết quả.

#### 4.2. Độ phức tạp và tính ổn định

> **Cần bổ sung:** Phân tích chi phí xây dựng ma trận và giải bài toán bình phương tối thiểu theo số quan sát $n$ và số điều hòa $K$. Độ phức tạp phụ thuộc phương pháp giải và kích thước ma trận (ma trận $A$ có kích thước $n \times (2K+1)$).

#### 4.3. Ưu điểm và hạn chế

**Ưu điểm cần thảo luận:**

- Mô hình có cấu trúc rõ ràng.
- Tuyến tính theo các hệ số khi chu kỳ đã biết.
- Có thể giải bằng đại số tuyến tính.
- Phù hợp với dữ liệu có tính chu kỳ.

**Hạn chế cần thảo luận:**

- Phụ thuộc vào chu kỳ và số điều hòa.
- Có thể không phù hợp với dữ liệu không tuần hoàn hoặc có thay đổi cấu trúc.
- Tăng số điều hòa có thể dẫn đến quá khớp.
- Chất lượng xấp xỉ phụ thuộc vào dữ liệu và cách lấy mẫu.

---

## PHẦN II. THỰC NGHIỆM BẰNG PYTHON

### 5. Mục tiêu và thiết kế thực nghiệm

#### 5.1. Câu hỏi thực nghiệm

- TLS xấp xỉ dữ liệu tuần hoàn tốt đến mức nào?
- Kết quả thay đổi ra sao khi số điều hòa $K$ tăng?
- Mô hình phản ứng thế nào khi dữ liệu có nhiễu?
- Kết quả có thay đổi khi chu kỳ giả định không chính xác?
- Mức độ phù hợp của TLS khác nhau như thế nào giữa dữ liệu tổng hợp và dữ liệu thực tế?

> **Cần bổ sung:** Chọn một số câu hỏi phù hợp với phạm vi và thời lượng tiểu luận.

#### 5.2. Môi trường và thư viện

- Python: [phiên bản]
- NumPy: xây dựng ma trận và giải bình phương tối thiểu.
- Matplotlib: trực quan hóa dữ liệu và đường cong xấp xỉ.
- Pandas (nếu cần): đọc và xử lý dữ liệu dạng bảng.
- SciPy (nếu cần): các công cụ phân tích tín hiệu hoặc hỗ trợ ước lượng chu kỳ.

> **Cần bổ sung:** Ghi phiên bản thư viện thực tế sau khi chạy chương trình.

### 6. Dữ liệu thực nghiệm

#### 6.1. Thực nghiệm A — Dữ liệu tổng hợp không nhiễu

Tạo dữ liệu từ một hàm lượng giác đã biết, ví dụ:

$$
y(x) = 2 + 1.5\cos\left(\frac{2\pi x}{T}\right) + 0.5\sin\left(\frac{4\pi x}{T}\right).
$$

- Khoảng giá trị của $x$: [...]
- Số điểm dữ liệu $n$: [...]
- Chu kỳ $T$: [...]
- Mục đích: kiểm tra xem TLS có khôi phục được các hệ số và tái tạo hàm đã cho hay không.

#### 6.2. Thực nghiệm B — Dữ liệu tổng hợp có nhiễu

Tạo dữ liệu quan sát:

$$
y_i = f(x_i) + \varepsilon_i,
$$

trong đó $\varepsilon_i$ là nhiễu được sinh theo một phân phối đã xác định, chẳng hạn nhiễu Gaussian có trung bình bằng 0 và độ lệch chuẩn $\sigma$.

- Mức nhiễu $\sigma$: [...]
- Hạt giống ngẫu nhiên (`random seed`): [...]
- Mục đích: đánh giá ảnh hưởng của nhiễu đến hệ số và độ chính xác xấp xỉ.

#### 6.3. Thực nghiệm C — Dữ liệu thực tế (tùy chọn)

> **Cần bổ sung:** Nếu sử dụng dữ liệu thực tế, mô tả nguồn dữ liệu, biến được phân tích, đơn vị đo, khoảng thời gian, tần suất lấy mẫu, cách xử lý giá trị thiếu và lý do dữ liệu có thể mang tính chu kỳ. Nếu không sử dụng dữ liệu thực tế, có thể bỏ mục này.

#### 6.4. Tiền xử lý và chia dữ liệu

> **Cần bổ sung:** Mô tả cách sắp xếp, làm sạch, xử lý giá trị thiếu và chuẩn hóa nếu cần. Nếu đánh giá khả năng dự báo, hãy phân chia theo thời gian thay vì xáo trộn ngẫu nhiên để tránh sử dụng thông tin tương lai. Phân biệt rõ xấp xỉ trên tập dữ liệu đã khớp với đánh giá ngoài mẫu.

### 7. Cài đặt TLS bằng Python

#### 7.1. Xây dựng ma trận thiết kế

> **Cần bổ sung:** Viết hàm nhận `x`, `T` và `K`, sau đó tạo các cột `1`, `cos(k*omega*x)` và `sin(k*omega*x)` với `k = 1, ..., K`. Kiểm tra kích thước ma trận và thứ tự các hệ số.

#### 7.2. Ước lượng hệ số

> **Cần bổ sung:** Sử dụng `numpy.linalg.lstsq(A, y, rcond=None)` hoặc một bộ giải bình phương tối thiểu phù hợp. Trình bày cách lấy vector hệ số và kiểm tra kết quả.

#### 7.3. Dự đoán và tính phần dư

> **Cần bổ sung:** Tính giá trị dự đoán tại các điểm quan sát và trên một lưới `x` mịn để vẽ đường cong. Tính phần dư `r_i = y_i - y_hat_i`.

#### 7.4. Mã nguồn minh họa

```python
import warnings

import numpy as np
import matplotlib.pyplot as plt


def true_function(x, period):
    """Underlying noise-free function used to generate synthetic data."""
    return (
        2.0
        + 1.5 * np.cos(2 * np.pi * x / period)
        + 0.5 * np.sin(4 * np.pi * x / period)
    )


def fit_tls(x, y, period, K):
    """Fit TLS for a known positive period.

    Returns coefficients ordered as [a0, a1, b1, ..., aK, bK].
    """
    x = np.asarray(x, dtype=float).ravel()
    y = np.asarray(y, dtype=float).ravel()

    if x.size != y.size or x.size == 0:
        raise ValueError("x and y must have the same non-zero length.")
    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
        raise ValueError("x and y must contain finite values only.")
    if not np.isfinite(period) or period <= 0:
        raise ValueError("period must be positive and finite.")
    if not isinstance(K, (int, np.integer)) or K < 0:
        raise ValueError("K must be a non-negative integer.")

    omega = 2 * np.pi / period
    columns = [np.ones_like(x)]
    for k in range(1, K + 1):
        columns.extend([
            np.cos(k * omega * x),
            np.sin(k * omega * x),
        ])

    A = np.column_stack(columns)
    beta, _, rank, _ = np.linalg.lstsq(A, y, rcond=None)

    if rank < A.shape[1]:
        warnings.warn("Design matrix is rank-deficient.")

    return beta


def predict_tls(x, beta, period):
    """Predict values using coefficients [a0, a1, b1, ..., aK, bK]."""
    x = np.asarray(x, dtype=float)
    beta = np.asarray(beta, dtype=float)
    omega = 2 * np.pi / period

    y_hat = np.full_like(x, beta[0], dtype=float)
    K = (len(beta) - 1) // 2
    for k in range(1, K + 1):
        a_k = beta[2 * k - 1]
        b_k = beta[2 * k]
        y_hat += (
            a_k * np.cos(k * omega * x)
            + b_k * np.sin(k * omega * x)
        )
    return y_hat


# Synthetic noisy data
rng = np.random.default_rng(42)
T = 24.0
x = np.linspace(0, 72, 180, endpoint=False)

y_true = true_function(x, T)
y_observed = y_true + rng.normal(0.0, 0.25, size=x.size)

K = 3
beta = fit_tls(x, y_observed, period=T, K=K)
y_fitted = predict_tls(x, beta, period=T)

rmse = np.sqrt(np.mean((y_observed - y_fitted) ** 2))
mae = np.mean(np.abs(y_observed - y_fitted))
print("Estimated coefficients:", beta)
print(f"Training RMSE: {rmse:.4f}")
print(f"Training MAE:  {mae:.4f}")

x_plot = np.linspace(x.min(), x.max(), 1000)
y_plot = predict_tls(x_plot, beta, period=T)
y_underlying = true_function(x_plot, T)

plt.figure(figsize=(10, 5))
plt.scatter(x, y_observed, s=18, label="Observed data")
plt.plot(x_plot, y_plot, label="TLS fit")
plt.plot(x_plot, y_underlying, "--", label="Underlying function")
plt.xlabel("x")
plt.ylabel("y")
plt.title(f"Trigonometric Least Squares (T={T}, K={K})")
plt.legend()
plt.tight_layout()
plt.show()
```

> **Cần bổ sung:** Chạy và kiểm tra mã nguồn trước khi đưa kết quả vào báo cáo. RMSE và MAE ở ví dụ này được tính trên chính dữ liệu dùng để khớp, vì vậy không được diễn giải chúng như sai số dự báo ngoài mẫu.

### 8. Chỉ số đánh giá

#### 8.1. Root Mean Squared Error (RMSE)

$$
\mathrm{RMSE} = \sqrt{\frac{1}{n} \sum_{i=1}^{n} (y_i - \hat y_i)^2}.
$$

RMSE nhấn mạnh hơn các sai số có độ lớn lớn do bình phương sai số trước khi lấy trung bình.

#### 8.2. Mean Absolute Error (MAE)

$$
\mathrm{MAE} = \frac{1}{n} \sum_{i=1}^{n} \left| y_i - \hat y_i \right|.
$$

MAE biểu diễn độ lớn sai số tuyệt đối trung bình và có cùng đơn vị với biến $y$.

#### 8.3. Sai số trên dữ liệu tổng hợp

> **Cần bổ sung:** Nếu biết hàm sinh dữ liệu, so sánh kết quả dự đoán với hàm nền không nhiễu. Có thể tính RMSE đối với $f(x_i)$ hoặc trên một lưới kiểm tra độc lập.

#### 8.4. Giới hạn của các chỉ số

> **Cần bổ sung:** Giải thích rằng sai số trên tập dữ liệu dùng để khớp thường giảm khi tăng số tham số, nhưng điều đó không đảm bảo khả năng tổng quát hóa tốt hơn. Khi so sánh các giá trị $K$, nên xem xét tập kiểm tra hoặc thiết kế đánh giá phù hợp.

### 9. Kết quả và thảo luận

#### 9.1. Kết quả trực quan

**Hình 1.** Dữ liệu tổng hợp và đường cong TLS.

> **Cần bổ sung:** Chèn hình và mô tả mức độ phù hợp của đường cong với dữ liệu.

**Hình 2.** Phần dư theo biến đầu vào.

> **Cần bổ sung:** Chèn biểu đồ phần dư và xem xét liệu còn cấu trúc chu kỳ rõ rệt hay không.

#### 9.2. Khảo sát số điều hòa $K$

> **Cần bổ sung:** Thử một tập giá trị $K$ phù hợp với kích thước dữ liệu, ví dụ $K \in \{1, 2, 3, 5, 8\}$ nếu ma trận vẫn đủ hạng. Lập bảng RMSE/MAE trên tập dữ liệu dùng để khớp và, nếu có, tập kiểm tra. Giải thích xu hướng quan sát được.

**Bảng 1.** Sai số theo số điều hòa $K$.

| Số điều hòa $K$ | RMSE (khớp) | MAE (khớp) | RMSE (kiểm tra) | Nhận xét |
|:---:|:---:|:---:|:---:|:---|
| … | … | … | … | … |
| … | … | … | … | … |
| … | … | … | … | … |

#### 9.3. Khảo sát mức nhiễu (tùy chọn)

> **Cần bổ sung:** Thay đổi độ lệch chuẩn của nhiễu và so sánh sai số, độ ổn định của các hệ số ước lượng.

#### 9.4. Khảo sát sai lệch chu kỳ (tùy chọn)

> **Cần bổ sung:** So sánh kết quả khi dùng chu kỳ đúng và một số chu kỳ giả định khác. Giữ các yếu tố khác cố định để phân tích ảnh hưởng của lựa chọn chu kỳ.

#### 9.5. Thảo luận tổng hợp

> **Cần bổ sung:** Trả lời trực tiếp từng câu hỏi thực nghiệm. Phân biệt kết quả quan sát được với giải thích giả thuyết; không đưa ra kết luận vượt quá phạm vi dữ liệu.

### 10. Kết luận

> **Cần bổ sung:** Tóm tắt nguyên lý TLS, phương pháp cài đặt, kết quả thực nghiệm nổi bật và các điều kiện mà phương pháp hoạt động tốt. Nêu hạn chế của thực nghiệm và hướng mở rộng, chẳng hạn ước lượng chu kỳ, so sánh với hồi quy đa thức hoặc Fourier, hay kiểm tra trên dữ liệu thực tế.

---

## Tài liệu tham khảo

> **Cần bổ sung:** Sách Giải tích số (Numerical Analysis), tài liệu về bình phương tối thiểu và xấp xỉ Fourier/lượng giác, cùng tài liệu chính thức của NumPy/Matplotlib nếu được sử dụng. Thống nhất một kiểu trích dẫn theo yêu cầu của giảng viên.

1. [Tác giả, năm, tên sách/bài báo, nhà xuất bản/tạp chí, thông tin xuất bản.]
2. [Tài liệu về bình phương tối thiểu và giải hệ tuyến tính.]
3. [Tài liệu về xấp xỉ lượng giác hoặc chuỗi Fourier.]
4. [Tài liệu chính thức của thư viện được sử dụng.]

## Phụ lục A. Mã nguồn

> **Cần bổ sung:** Đưa mã nguồn hoàn chỉnh hoặc liên kết kho mã nếu giảng viên cho phép.

## Phụ lục B. Bảng kết quả bổ sung

> **Cần bổ sung:** Đưa các bảng số liệu, cấu hình thực nghiệm hoặc hình ảnh bổ sung.

---

## Checklist trước khi nộp

- [ ] Đã định nghĩa đầy đủ ký hiệu và giả thiết.
- [ ] Đã phân biệt rõ chu kỳ $T$ và số điều hòa $K$.
- [ ] Đã giải thích cách xây dựng ma trận thiết kế.
- [ ] Đã dùng bộ giải bình phương tối thiểu phù hợp và kiểm tra hạng ma trận.
- [ ] Đã chạy lại mã nguồn và kiểm tra các kết quả được báo cáo.
- [ ] Đã ghi rõ dữ liệu, tham số, phiên bản thư viện và cách tái lập thực nghiệm.
- [ ] Đã phân biệt sai số trên dữ liệu khớp với sai số trên dữ liệu kiểm tra.
- [ ] Đã đánh số hình, bảng, phương trình và trích dẫn tài liệu nhất quán.
- [ ] Đã ghi nhận các hạn chế thay vì chỉ trình bày ưu điểm của phương pháp.
