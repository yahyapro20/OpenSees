# راهنمای Build OpenSees — Intel MKL + MUMPS

**Build Guide v3.0 — Native Linux**

راهنمای کامل Build کردن OpenSees با Intel MKL و MUMPS Solver

نسخه High-Performance، بهینه‌شده برای پردازنده‌های چند‌هسته‌ای مدرن با AVX2 — مستقیماً روی لینوکس، بدون واسطه WSL یا Windows.

| مشخصات | مقدار |
|---|---|
| 🖥️ سیستم مرجع تست | Intel Core i9-13900K 32GB DDR5 — Ubuntu 22.04 LTS |
| 🧵 هسته‌های منطقی | 24 Threads |
| ⚙️ Compiler | Intel ifx + GCC |
| 🧮 Math Library | Intel MKL 2026 |
| 🔗 Solver | MUMPS 5.5.1 |

> ⚠️ AVX-512 روی نسل i9-13900K (Raptor Lake) وجود ندارد. Intel این قابلیت را در این نسل غیرفعال کرده — فقط AVX2 واقعی است و هیچ compiler flag ای آن را فعال نمی‌کند. اگر CPU دیگری داری، قبل از استفاده از `-mavx512f` حتماً پشتیبانی واقعی را با `lscpu` چک کن.

---

## فهرست مراحل

**بخش اول: Build Guide**

1. [آماده‌سازی سیستم](#1-آماده‌سازی-سیستم-و-نصب-dependencies)
2. [نصب Intel oneAPI](#2-نصب-intel-oneapi--mkl--mpi--fortran)
3. [Build کردن MUMPS](#3-build-کردن-mumps)
4. [Clone و Build OpenSees](#4-clone-و-build-کردن-opensees--openseessp)
5. [تست و راه‌اندازی نهایی](#5-تست-و-راه‌اندازی-نهایی)
6. [خطاها و راه‌حل‌ها](#6-جدول-خطاها-و-راه‌حل‌ها)

**بخش دوم: مکانیزم تحلیل**

7. [اجزای Analysis](#7-اجزای-تحلیل-و-ترتیب-تعریف-آن‌ها)
8. [جریان تحلیل خطی](#8-جریان-اجرای-تحلیل-خطی-linear-static)
9. [جریان تحلیل غیرخطی](#9-جریان-اجرای-تحلیل-غیرخطی-newton-raphson)

---

## 1. آماده‌سازی سیستم و نصب Dependencies

⏱️ ۱۰ دقیقه | 🟢 آسان

این راهنما فرض می‌کند از قبل یک سیستم Ubuntu 22.04 (یا هر توزیع سازگار دیگر) در اختیار داری — چه به‌صورت native نصب‌شده، چه روی یک VM یا سرور — و به آن دسترسی sudo داری. نیازی به WSL یا هیچ ابزار مخصوص Windows نیست؛ تمام دستورات مستقیماً در ترمینال لینوکس اجرا می‌شوند.

> 💡 **قبل از شروع:** ترمینال یک پنجره‌ی متنی برای دادن دستور مستقیم به سیستم‌عامل است — همان کاری که در ادامه‌ی این راهنما مدام در آن دستور کپی/پیست می‌کنی. `sudo` یعنی «این دستور را با دسترسی مدیر سیستم اجرا کن» — چیزهایی مثل نصب نرم‌افزار جدید بدون آن ممکن نیست. هر جا `sudo` دیدی، سیستم رمز عبور کاربری‌ات را می‌پرسد.

### به‌روزرسانی سیستم

```bash
sudo apt update && sudo apt upgrade -y
```

```
# خروجی نمونه
Reading package lists... Done
Building dependency tree... Done
0 upgraded, 0 newly installed, 0 to remove.
```

> 💡 **apt یعنی چه؟** `apt` مدیر پکیج (package manager) اوبونتو است — انبار مرکزی نرم‌افزارهایی که تیم اوبونتو از قبل کامپایل و تست کرده. `update` فهرست این انبار را تازه می‌کند، `upgrade` نرم‌افزارهای نصب‌شده را به آخرین نسخه می‌رساند. تا آخر این راهنما، هر جا `sudo apt install X` ببینی یعنی «این نرم‌افزار آماده را از انبار رسمی نصب کن» — نه اینکه خودمان از صفر بسازیمش.

### نصب Build Tools و Dependencies

> ✅ همه پکیج‌های زیر را در همین مرحله نصب کن. این لیست کامل است و شامل تمام پکیج‌هایی است که اگر نباشند، CMake در مراحل بعد خطا می‌دهد.

```bash
sudo apt install -y \
     build-essential \
     gcc g++ gfortran \
     git \
     ninja-build \
     pkg-config \
     wget curl \
     python3 python3-pip \
     liblapack-dev libblas-dev \
     libopenmpi-dev openmpi-bin \
     libhdf5-dev \
     tcl-dev \
     tk-dev \
     libtk8.6 \
     libeigen3-dev \
     libssl-dev \
     numactl \
     hwloc
```

**دلیل هر پکیج حیاتی:**

| پکیج | دلیل | خطا در صورت نبود |
|---|---|---|
| `libhdf5-dev` | HDF5 recorder | `MPCORecorder.cpp: hdf5_hl.h not found` |
| `tk-dev` + `libtk8.6` | رابط گرافیکی Tcl | `Could NOT find TCLTK` |
| `libeigen3-dev` | کتابخانه ماتریس | `Could not find Eigen3` |
| `libopenmpi-dev` | MPI برای اجرای موازی | `MPI not found` |

> 💡 **Compile کردن یعنی چه؟** OpenSees به زبان‌های C و C++ نوشته شده — کدی که فقط انسان و کامپایلر آن را می‌فهمند، نه CPU. «Compile کردن» یعنی این کد متنی را با یک برنامه به‌نام کامپایلر (اینجا: `gcc`/`g++` برای C/C++ و `gfortran` بعداً `ifx` برای Fortran) به دستورات ماشین تبدیل کنیم. بعد از آن، برنامه‌ای به‌نام لینکر این قطعات را همراه با کتابخانه‌های آماده (مثل همین پکیج‌های `lib...-dev` بالا) به یک فایل اجرایی نهایی می‌چسباند. پسوند `-dev` یعنی «نسخه‌ی توسعه» یک کتابخانه — همان فایل‌هایی که کامپایلر و لینکر برای ساختن یک برنامه‌ی جدید از روی آن کتابخانه لازم دارند (نسخه‌ی معمولی بدون `-dev` فقط برای اجرای برنامه‌های آماده کافی است، نه ساختن برنامه‌ی جدید).

### نصب CMake جدید (نسخه apt کافی نیست)

> ⚠️ نسخه CMake در apt ریپازیتوری Ubuntu 22.04 برابر 3.22 است، اما OpenSees به 3.24 به بالا نیاز دارد. باید از ریپازیتوری رسمی Kitware نصب شود.

```bash
# حذف نسخه قدیمی
sudo apt remove cmake -y

# اضافه کردن Kitware repo
wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc 2>/dev/null | \
    gpg --dearmor - | \
    sudo tee /usr/share/keyrings/kitware-archive-keyring.gpg >/dev/null

echo 'deb [signed-by=/usr/share/keyrings/kitware-archive-keyring.gpg] https://apt.kitware.com/ubuntu/ jammy main' | \
    sudo tee /etc/apt/sources.list.d/kitware.list >/dev/null

sudo apt update
sudo apt install cmake -y

# تأیید — باید 3.31+ باشد
cmake --version
```

> 💡 **CMake چه‌کاری انجام می‌دهد؟** CMake خودش کامپایلر نیست. پروژه‌ی OpenSees صدها فایل کد دارد؛ CMake با خواندن دستورالعمل‌های پروژه (فایل‌های `CMakeLists.txt`) مشخص می‌کند کدام فایل‌ها باید به چه ترتیبی، با چه فلگ‌هایی، و با چه کتابخانه‌هایی کامپایل و لینک شوند — و این دستورالعمل نهایی را برای ابزار build واقعی (اینجا معمولاً `make` یا `ninja`) تولید می‌کند. به همین دلیل بهش «build system generator» می‌گویند، نه خود سیستم build.

### نصب Conan

```bash
pip3 install conan
echo 'export PATH="$HOME/.local/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc
conan --version
```

> 💡 **چرا هم apt هم Conan؟** `apt` کتابخانه‌های سطح سیستم‌عامل را مدیریت می‌کند (چیزهایی که خیلی از برنامه‌ها به آن نیاز دارند). Conan یک مدیر پکیج مخصوص پروژه‌های C++ است — بعضی کتابخانه‌های خاصی که OpenSees نیاز دارد نسخه‌ی apt ندارند یا نسخه‌شان قدیمی است؛ CMake در مرحله‌ی بعد از Conan می‌خواهد آن‌ها را دقیقاً با نسخه‌ی درست دانلود و آماده کند.

---

## 2. نصب Intel oneAPI — MKL + MPI + Fortran

⏱️ ۱۵ دقیقه | 🟡 متوسط

این مرحله حدود ۳ گیگابایت دانلود دارد و شامل Intel MKL (کتابخانه جبر خطی)، Intel MPI و کامپایلر Fortran (`ifx`) می‌شود.

> 💡 **چرا اصلاً MKL لازم است؟** هسته‌ی محاسباتی OpenSees در هر گام تحلیل، حل یک دستگاه معادلات خطی بزرگ است — عملیات جبر خطی (ضرب و تجزیه‌ی ماتریس). MKL (Math Kernel Library) نسخه‌ای از این عملیات است که Intel آن را مخصوص CPUهای خودش به‌شدت بهینه کرده (استفاده از دستورات SIMD مثل AVX2). بدون آن، OpenSees از یک پیاده‌سازی عمومی و کند دهه‌ی ۹۰ (Netlib LAPACK) استفاده می‌کند که ده‌ها برابر کندتر است — برنامه درست کار می‌کند، فقط تحلیل‌های بزرگ ساعت‌ها بیشتر طول می‌کشد.

```bash
# اضافه کردن Intel repo
wget -O- https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB | \
    gpg --dearmor | \
    sudo tee /usr/share/keyrings/oneapi-archive-keyring.gpg > /dev/null

echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" | \
    sudo tee /etc/apt/sources.list.d/oneAPI.list

sudo apt update

# نصب پکیج‌های لازم
sudo apt install -y \
    intel-oneapi-mkl \
    intel-oneapi-mkl-devel \
    intel-oneapi-compiler-fortran-2025.1 \
    intel-oneapi-mpi \
    intel-oneapi-mpi-devel
```

### تنظیم PATH و Library Path

> ⚠️ بعد از نصب، `ifx` در PATH نیست و `libifcoremt.so.5` پیدا نمی‌شود. هر دو باید دستی تنظیم شوند وگرنه OpenSees اجرا نمی‌شود.

> 💡 **دو مفهوم پشت این خطا:** PATH یک متغیر محیطی (environment variable) است — لیستی از پوشه‌ها که ترمینال هنگام تایپ یک دستور (مثل `ifx`) در آن‌ها دنبال فایل اجرایی می‌گردد. چون Intel آن را در یک پوشه‌ی غیر پیش‌فرض نصب کرده، تا این پوشه به PATH اضافه نشود سیستم نمی‌داند `ifx` کجاست. مورد دوم فرق دیگری دارد: `libifcoremt.so.5` یک کتابخانه‌ی پویا (dynamic library) است — کدی که به‌جای اینکه داخل هر برنامه کپی شود، فقط هنگام اجرا از دیسک بارگذاری می‌شود (برخلاف کتابخانه‌ی استاتیک که بعداً می‌بینیم). سیستم برای پیدا کردن این‌جور فایل‌ها یک مسیر جدا به‌نام `ld.so.conf.d` دارد که با `ldconfig` به‌روزرسانی می‌شود — به همین دلیل این خطا با اضافه‌کردن به PATH ساده حل نمی‌شود.

```bash
# اضافه کردن ifx به PATH
echo 'export PATH="/opt/intel/oneapi/compiler/2025.1/bin:$PATH"' >> ~/.bashrc

# فعال‌سازی خودکار oneAPI در هر session
echo 'source /opt/intel/oneapi/setvars.sh --force > /dev/null 2>&1' >> ~/.bashrc

# تنظیم دائمی dynamic library path
# (بدون این: libifcoremt.so.5: cannot open shared object file)
echo '/opt/intel/oneapi/compiler/2025.1/lib' | sudo tee /etc/ld.so.conf.d/intel-oneapi.conf
sudo ldconfig

source ~/.bashrc
```

### تأیید نهایی محیط

```bash
echo "=== Build Environment ==="
echo "GCC:    $(gcc --version | head -1)"
echo "GFortran: $(gfortran --version | head -1)"
echo "ifx:    $(ifx --version | head -1)"
echo "CMake:  $(cmake --version | head -1)"
echo "MPI:    $(mpirun --version 2>&1 | head -1)"
echo "MKL:    $MKLROOT"
echo "Conan:  $(conan --version)"
```

```
# خروجی تأیید‌شده
GCC:    gcc (Ubuntu 11.4.0-1ubuntu1~22.04.3) 11.4.0
GFortran: GNU Fortran (Ubuntu 11.4.0-1ubuntu1~22.04.3) 11.4.0
ifx:    ifx (IFX) 2025.1.1 20250418
CMake:  cmake version 4.3.2
MPI:    Intel(R) MPI Library for Linux* OS, Version 2021.18.0
MKL:    /opt/intel/oneapi/mkl/2026.0
Conan:  Conan version 2.28.1
```

---

## 3. Build کردن MUMPS

⏱️ ۱۰ دقیقه | 🟡 متوسط

> 💡 **MUMPS چیست و چرا جدا build می‌شود؟** MUMPS (MUltifrontal Massively Parallel Sparse direct Solver) یک کتابخانه‌ی حل معادلات خطی موازی و توزیع‌شده است — یعنی می‌تواند حل یک دستگاه معادلات بزرگ را بین چند پردازنده/سرور تقسیم کند (برخلاف MKL که روی یک ماشین چندهسته‌ای کار می‌کند). چون MUMPS پکیج آماده در apt ندارد و OpenSeesSP مستقیماً به فایل‌های کامپایل‌شده‌ی آن نیاز دارد، باید قبل از خود OpenSees، آن را جدا از سورس کامپایل کنیم.

> 🚨 آدرس رسمی `http://mumps.enseeiht.fr/MUMPS_5.5.1.tar.gz` از دسترس خارج است. از GitHub mirror استفاده می‌کنیم.

```bash
mkdir -p ~/opensees-build
cd ~/opensees-build

# clone از GitHub mirror (کار می‌کند)
git clone https://github.com/scivision/mumps.git mumps
cd mumps

mkdir build && cd build
source /opt/intel/oneapi/setvars.sh --force 2>/dev/null

cmake .. \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_Fortran_COMPILER=ifx \
  -DCMAKE_C_COMPILER=gcc \
  -DCMAKE_CXX_COMPILER=g++ \
  -DMUMPS_parallel=ON \
  -DMUMPS_openmp=ON \
  -DCMAKE_INSTALL_PREFIX=~/opensees-build/mumps/install

cmake --build . --parallel 24
cmake --install .
```

> ℹ️ ممکن است در مرحله test خطای `-fiopenmp` بدهد — این مشکلی نیست. کتابخانه‌های اصلی MUMPS نصب می‌شوند.

### تأیید MUMPS

```bash
ls ~/opensees-build/mumps/install/lib/
```

```
# خروجی موفق
libdmumps.a  libmumps_common.a  libpord.a  libsmumps.a
```

> 💡 **پسوند `.a` یعنی چه؟** فایل‌های `.a` که همین الان ساخته شدند کتابخانه‌ی استاتیک (static library) هستند — بسته‌ای از کدهای از‌قبل‌کامپایل‌شده که، برخلاف کتابخانه‌ی پویا (`.so`) که در مرحله‌ی قبل دیدیم، هنگام لینک کردن مستقیماً داخل فایل اجرایی نهایی کپی می‌شوند. یعنی OpenSeesSP که بعداً می‌سازیم، کد MUMPS را همراه خودش حمل می‌کند و در زمان اجرا به فایل جدایی از MUMPS نیاز ندارد.

---

## 4. Clone و Build کردن OpenSees + OpenSeesSP

⏱️ ۱۰ دقیقه | 🟡 متوسط

### Clone کردن OpenSees

```bash
cd ~/opensees-build
git clone https://github.com/OpenSees/OpenSees.git
cd OpenSees
```

> 💡 **git clone یعنی چه؟** Git ابزار مدیریت نسخه‌ی کد است و GitHub سروری است که سورس‌کد OpenSees روی آن نگه‌داری می‌شود. دستور `git clone` یک کپی کامل از تمام فایل‌های سورس (و تاریخچه‌ی تغییراتش) را از GitHub روی سیستم تو دانلود می‌کند — این کد متنی هنوز کامپایل نشده؛ همان چیزی است که در مراحل بعد قرار است با CMake و کامپایلر به یک برنامه‌ی اجرایی تبدیلش کنیم.

### CMake Configure

```bash
mkdir build && cd build
source /opt/intel/oneapi/setvars.sh --force 2>/dev/null

cmake .. \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER=gcc \
  -DCMAKE_CXX_COMPILER=g++ \
  -DCMAKE_Fortran_COMPILER=ifx \
  -DCMAKE_C_FLAGS="-O3 -march=native -mavx2 -flto" \
  -DCMAKE_CXX_FLAGS="-O3 -march=native -mavx2 -flto" \
  -DCMAKE_Fortran_FLAGS="-O3 -march=native -mavx2 -flto" \
  -DBLA_VENDOR=Intel10_64lp \
  -DMKL_LINK=static \
  -DMKL_INTERFACE_FULL=intel_lp64 \
  -DMUMPS_DIR="$HOME/opensees-build/mumps/install" \
  -DCMAKE_INSTALL_PREFIX="$HOME/opensees-build/install"

echo "=== CMake exit code: $? ==="
```

> 🚨 عمداً از `-ffast-math` استفاده نشده. این flag فرض می‌کند NaN/Inf وجود ندارند و ترتیب عملیات اعشاری را برای سرعت بازآرایی می‌کند — یعنی نتیجه دیگر دقیقاً IEEE‑754 نیست. در تحلیل غیرخطی که به iteration‌های Newton-Raphson و چک‌های واگرایی (NaN/Inf) وابسته است، این می‌تواند مسیر همگرایی یا حتی نتیجه نهایی را بی‌سروصدا تغییر دهد. سرعتی که از دست می‌رود ناچیز است چون بخش سنگین محاسبات (فاکتورگیری ماتریس) داخل MKL/MUMPS انجام می‌شود، نه در کد OpenSees.

**توضیح flags مهم:**

| Flag | توضیح |
|---|---|
| `-O3` | حداکثر سطح بهینه‌سازی compiler — کاملاً IEEE-safe |
| `-march=native` | فعال‌سازی تمام قابلیت‌های CPU میزبان (مثل AVX2، FMA، BMI2) |
| `-mavx2` | فعال‌سازی صریح AVX2 (SIMD ۲۵۶‑بیتی) |
| `-flto` | Link-Time Optimization — inlining بین فایل‌ها، بدون تغییر نتیجه محاسبات |
| `Intel10_64lp` | MKL با LP64 interface استاندارد |
| `static` | MKL به‌صورت static لینک می‌شود |

```
# خروجی موفق CMake
-- Found BLAS: /opt/intel/oneapi/mkl/2026.0/lib/libmkl_intel_lp64.so ...
-- Found LAPACK: ...
-- Found HDF5: ... (found version "1.10.7")
-- Found TCLTK: /usr/lib/x86_64-linux-gnu/libtcl.so
-- MKL was found.
-- MPI was found.
-- Configuring done
-- Generating done
```

> 💡 **target یعنی چه؟** یک پروژه‌ی CMake می‌تواند از یک سورس‌کد، چند فایل اجرایی متفاوت بسازد — هرکدام را یک target می‌گویند. مرحله‌ی قبل (Configure) فقط دستورالعمل build را آماده کرد؛ الان با `cmake --build . --target <نام>` مشخص می‌کنیم دقیقاً کدام target را می‌خواهیم بسازیم. OpenSees دو target مجزا دارد: `OpenSees` (نسخه‌ی سریال، تک‌پردازنده) و `OpenSeesSP` (نسخه‌ای که با MPI روی چند پردازنده/سرور به‌صورت موازی اجرا می‌شود) — و باید هرکدام را جدا build کرد.

### Build OpenSees

```bash
cmake --build . --target OpenSees --parallel 24
echo "Exit: $?"
```

> ℹ️ با ۲۴ هسته موازی، build حدود ۲ تا ۳ دقیقه طول می‌کشد. اگر تعداد هسته‌های تو کمتر است، زمان build به همان نسبت بیشتر می‌شود — کاملاً طبیعی است.

### Build OpenSeesSP (نسخه Parallel)

> 🚨 **این بخش حیاتی است.** target بالا (`OpenSees`) فقط باینری سریال را می‌سازد. اگر از OpenSeesSP استفاده می‌کنی، باید این target جدا را هم build کنی — بدون آن، تمام تنظیمات MUMPS/MPI که در مراحل قبل انجام دادیم عملاً بلااستفاده می‌ماند.

> 💡 **فرق OpenSees و OpenSeesSP:** نسخه‌ی سریال (`OpenSees`) کل مدل و تحلیل را روی یک پردازنده اجرا می‌کند. نسخه‌ی موازی (`OpenSeesSP`) با استفاده از MPI (Message Passing Interface) مدل را بین چند «فرآیند» (process) — که می‌توانند حتی روی چند کامپیوتر جدا باشند — تقسیم می‌کند؛ هر فرآیند سهم خودش را حل می‌کند و از طریق MPI با بقیه پیام رد و بدل می‌کند. این یعنی OpenSeesSP یک برنامه‌ی جداست، نه یک flag اضافه روی همان برنامه‌ی سریال — به همین دلیل target مجزا دارد.

```bash
cmake --build . --target OpenSeesSP --parallel 24
echo "Exit: $?"
ls -la ./OpenSeesSP
```

### تأیید: آیا واقعاً MKL لینک شده یا fallback به Netlib قدیمی؟

> ⚠️ اگر `find_package(LAPACK)` به هر دلیلی MKL را پیدا نکند، CMake بدون هیچ خطایی به BLAS/LAPACK رفرنس دهه ۹۰ (نسخه apt اوبونتو) سقوط می‌کند. build موفق می‌شود، OpenSeesSP هم اجرا می‌شود و نتیجه هم درست است — فقط ده‌ها برابر کندتر، بدون هیچ هشداری. دو راه برای اطمینان:
>
> 1. لاگ CMake configure (چند بخش بالاتر) — خط `Found BLAS` باید مسیر `/opt/intel/oneapi/mkl/...` را نشان دهد. اگر به‌جایش `/usr/lib/x86_64-linux-gnu/libblas.so.3` یا مشابه دیدی، یعنی fallback شده.
> 2. چک باینری نهایی — چون MKL را static لینک کرده‌ایم، `ldd` کمکی نمی‌کند؛ باید داخل باینری دنبال رشته نسخه MKL بگردیم.

```bash
strings ./OpenSeesSP | grep -i "Math Kernel Library"
```

```
# اگر MKL واقعاً لینک شده باشد
Intel(R) Math Kernel Library Version 2026.0.0 Product Build ...
```

> 🚨 اگر این دستور هیچ خروجی‌ای نداد، یعنی به BLAS/LAPACK رفرنس fallback شده — برگرد و مطمئن شو `source /opt/intel/oneapi/setvars.sh` قبل از `cmake ..` اجرا شده و `MKLROOT` ست است، بعد `rm -rf build` و configure را از نو انجام بده.

---

## 5. تست و راه‌اندازی نهایی

⏱️ ۸ دقیقه | 🟢 آسان

### تست اجرا

```bash
cd ~/opensees-build/OpenSees/build
echo 'puts {OpenSees is working}; exit' | ./OpenSees
```

```
# خروجی موفق
   OpenSees -- Open System For Earthquake Engineering Simulation
        Pacific Earthquake Engineering Research Center
Version 3.8.0 64-Bit
(c) Copyright 1999-2016 The Regents of the University of California

OpenSees is working
```

### تست OpenSeesSP + MUMPS

```bash
echo 'system Mumps; puts {MUMPS solver OK}; exit' > /tmp/test_mumps.tcl
mpirun -np 2 ./OpenSeesSP /tmp/test_mumps.tcl
```

```
# خروجی موفق
MUMPS solver OK
```

> 💡 **mpirun -np 2 یعنی چه؟** `mpirun` راه‌انداز برنامه‌های MPI است. فلگ `-np 2` یعنی «این برنامه را در ۲ فرآیند موازی اجرا کن» — همان مفهوم فرآیندهای موازی که در مرحله‌ی قبل دیدیم. برای تحلیل واقعی، این عدد را معمولاً برابر یا کمتر از تعداد هسته‌های آزاد سیستم می‌گذاری؛ عددی بیشتر از تعداد هسته‌ها معمولاً باعث کند شدن به‌جای سریع شدن می‌شود.

> 🚨 **نکته حیاتی برای استفاده واقعی:** هرچقدر هم build درست باشد، اگر توی Tcl script تحلیل خودت دستور `system Mumps` را ننویسی (و به‌جایش solver پیش‌فرض بماند)، از حل موازی/توزیع‌شده استفاده نمی‌شود. این خط باید قبل از `analyze` و بعد از تعریف `constraints`/`numberer` در اسکریپت بیاید — مثلاً:

```tcl
constraints Transformation
numberer ParallelRCM
system Mumps          ;# این خط — بدونش solver پیش‌فرض جایگزین می‌شود
test NormDispIncr 1.0e-6 25
algorithm Newton
integrator LoadControl 0.01
analysis Static
```

### نصب دائمی و تنظیم PATH

```bash
mkdir -p ~/opensees-build/install/bin
cp ~/opensees-build/OpenSees/build/OpenSees ~/opensees-build/install/bin/
echo 'export PATH="$HOME/opensees-build/install/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc

# تست از هر جا
OpenSees --version
```

### چک‌لیست نهایی

- [x] سیستم Ubuntu 22.04 (یا معادل) به‌روز و آماده است
- [x] CMake نسخه 3.31 به بالا نصب شده
- [x] Intel MKL 2026.0 نصب و تنظیم شده
- [x] ifx (Intel Fortran) در PATH است
- [x] libifcoremt.so.5 در ldconfig ثبت شده
- [x] MUMPS build شده (libdmumps.a موجود است)
- [x] OpenSees با `-O3 -march=native -mavx2` build شده
- [x] OpenSees اجرا می‌شود و پیام موفق نشان می‌دهد
- [x] target جدای OpenSeesSP هم build شده (نه فقط OpenSees سریال)
- [x] با `strings` روی باینری تأیید شده که MKL واقعی لینک شده، نه fallback به Netlib
- [x] OpenSeesSP با `mpirun` و `system Mumps` تست شده و بدون خطا اجرا شده
- [x] اسکریپت‌های Tcl تحلیل واقعی شامل دستور `system Mumps` هستند

---

## 6. جدول خطاها و راه‌حل‌ها

| خطا | علت | راه‌حل |
|---|---|---|
| `hdf5_hl.h not found` | `libhdf5-dev` نصب نیست | `sudo apt install libhdf5-dev` |
| `Could NOT find TCLTK` | `tk-dev` نصب نیست | `sudo apt install tk-dev libtk8.6` |
| `Could not find Eigen3` | `libeigen3-dev` نصب نیست | `sudo apt install libeigen3-dev` |
| `MUMPS download failed` | آدرس رسمی از دسترس خارج است | `git clone https://github.com/scivision/mumps.git` |
| `conan: command not found` | مسیر pip در PATH نیست | `echo 'export PATH="$HOME/.local/bin:$PATH"' >> ~/.bashrc` |
| `ifx: command not found` | مسیر Intel compiler در PATH نیست | `export PATH="/opt/intel/oneapi/compiler/2025.1/bin:$PATH"` |
| `libifcoremt.so.5: cannot open` | مسیر runtime library تنظیم نشده | `echo '/opt/intel/oneapi/compiler/2025.1/lib' \| sudo tee /etc/ld.so.conf.d/intel-oneapi.conf && sudo ldconfig` |
| کاراکترهای عجیب هنگام paste (مثل `^[[200~`) | bracketed-paste در بعضی ترمینال‌ها کاراکتر اضافه می‌آورد | دستورات چندخطی را به‌جای paste مستقیم، تایپ کن یا از یک فایل اسکریپت اجرا کن |
| CMake نسخه 3.22 گزارش می‌شود | نسخه قدیمی apt در Ubuntu | نصب از ریپازیتوری رسمی Kitware (مرحله ۱) |
| CMake cache خراب / ناسازگار | پوشه build قدیمی باقی مانده | `rm -rf build && mkdir build && cd build` |
| `-fiopenmp` خطا در تست MUMPS | GCC این flag مخصوص Intel را نمی‌شناسد | نادیده بگیر — کتابخانه‌های اصلی همچنان نصب می‌شوند |

---

## بخش دوم — مفاهیم تحلیل

**OpenSees دقیقاً چطور یک تحلیل را حل می‌کند؟**

این بخش دیگر درباره‌ی build کردن نیست — درباره‌ی این است که وقتی اسکریپت Tcl دستور `analyze` را صدا می‌زند، پشت‌صحنه‌ی OpenSees دقیقاً چه اتفاقی، به چه ترتیبی، می‌افتد؛ هم برای تحلیل خطی و هم غیرخطی.

---

## 7. اجزای تحلیل و ترتیب تعریف آن‌ها

⏱️ ۱۰ دقیقه مطالعه | 📘 مفهومی

قبل از `analyze`، مدل (گره‌ها، المان‌ها، بارها) با دستوراتی مثل `model`، `node`، `element` و `pattern` در حافظه‌ای به‌نام Domain ساخته می‌شود. اما این‌که این مدل چطور حل شود، به ۶ جزء دیگر بستگی دارد که باید دقیقاً به همین ترتیب، پیش از `analyze`، در اسکریپت تعریف شوند:

| ترتیب | دستور | نقش | مثال رایج |
|---|---|---|---|
| ۱ | `constraints` | نحوه‌ی اعمال محدودیت‌های مرزی (تکیه‌گاه‌ها، هم‌بستگی گره‌ها) | `Transformation`, `Plain`, `Penalty` |
| ۲ | `numberer` | ترتیب شماره‌گذاری درجات آزادی، برای کارایی حل معادلات | `RCM`, `ParallelRCM` (SP) |
| ۳ | `system` | نوع دستگاه معادلات خطی و solverی که آن را حل می‌کند | `BandGeneral`, `SparseGeneral`, `Mumps` |
| ۴ | `test` | معیار همگرایی — فقط در تحلیل غیرخطی معنا دارد | `NormDispIncr`, `EnergyIncr` |
| ۵ | `algorithm` | روش حل معادلات در هر گام (تکراری یا یک‌باره) | `Linear`, `Newton`, `ModifiedNewton` |
| ۶ | `integrator` | نحوه‌ی پیش‌روی بار یا زمان بین گام‌ها | `LoadControl`, `Newmark` |
| ۷ | `analysis` | نوع کلی تحلیل — این دستور بقیه‌ی اجزا را که از قبل تعریف شده‌اند به هم متصل می‌کند | `Static`, `Transient` |

> 💡 **چرا ترتیب مهم است؟** دستور `analysis` در انتها، تمام اجزای بالا را که از قبل تعریف شده‌اند به هم متصل می‌کند. اگر یکی از این‌ها قبل از `analysis` تعریف نشده باشد، OpenSees مقدار پیش‌فرض داخلی خودش را جایگزین می‌کند — بدون هیچ خطایی. دقیقاً همان چیزی که در راهنمای build (مرحله ۵) درباره‌ی `system Mumps` دیدیم: نبودنش خطا نمی‌دهد، فقط بی‌صدا solver را عوض می‌کند.

```tcl
# ترتیب صحیح — دقیقاً همین توالی
constraints Transformation
numberer RCM
system BandGeneral
test NormDispIncr 1.0e-6 25
algorithm Newton
integrator LoadControl 0.01
analysis Static

# فقط الان که همه‌چیز تعریف شده، analyze معنا دارد
analyze 100
```

---

## 8. جریان اجرای تحلیل خطی (Linear Static)

⏱️ ۵ دقیقه مطالعه | 📘 مفهومی

وقتی `algorithm Linear` انتخاب شده، به‌ازای هر گامی که `analyze` اجرا می‌کند، دقیقاً همین توالی — یک‌بار، بدون هیچ تکراری — رخ می‌دهد:

1. **اعمال بار گام فعلی** — integrator مقدار بار این گام را به Domain اضافه می‌کند (مثلاً طبق LoadControl)
2. **تشکیل ماتریس سختی K** — یک‌بار، بر اساس رفتار خطی فرض‌شده‌ی مصالح و هندسه
3. **حل معادله — `K · Δu = ΔP`** — توسط solver مشخص‌شده در system (مثلاً BandGeneral یا Mumps)
4. **به‌روزرسانی جابجایی‌های گرهی** — Δu به وضعیت فعلی مدل اضافه می‌شود
5. ✅ **پایان گام** — بدون هیچ بررسی همگرایی. چون K ثابت فرض شده، جواب یک‌بار حل‌شدن دقیقاً همان جواب نهایی است. حتی اگر test هم تعریف شده باشد، در این مسیر اصلاً استفاده نمی‌شود.

> ℹ️ به همین دلیل تحلیل خطی همیشه سریع‌تر از غیرخطی است: به‌ازای هر گام فقط یک‌بار ماتریس تشکیل و حل می‌شود، نه چندین تکرار.

---

## 9. جریان اجرای تحلیل غیرخطی (Newton-Raphson)

⏱️ ۱۰ دقیقه مطالعه | 📘 مفهومی

وقتی مصالح یا هندسه رفتار غیرخطی دارند، ماتریس سختی K دیگر ثابت نیست — به وضعیت لحظه‌ای مدل بستگی دارد. به همین دلیل، برخلاف حالت خطی، هر گام بار باید با تکرار (iteration) حل شود تا به جواب واقعی برسد. این چرخه را Newton-Raphson می‌گویند، و اینجاست که `algorithm` و `test` نقش اصلی را ایفا می‌کنند:

1. **اعمال افزایش بار گام** — integrator (مثلاً `LoadControl 0.01`) بار جدید را اضافه می‌کند
2. **تشکیل ماتریس سختی مماسی (tangent)** — بر اساس وضعیت *فعلی* تنش/کرنش هر المان
3. **حل معادله برای Δu** — توسط solver مشخص‌شده در system
4. **به‌روزرسانی جابجایی‌ها و محاسبه‌ی residual** — نیروی نامتعادل باقی‌مانده = بار اعمال‌شده − مقاومت داخلی فعلی مدل
5. ❓ **test بررسی می‌کند:** آیا residual/جابجایی افزایشی از تلورانس تعریف‌شده کوچک‌تر شده؟ (مثلاً `NormDispIncr 1.0e-6 25`)
   - ⟲ اگر همگرا نشده و به سقف iteration (اینجا ۲۵) نرسیده‌ایم → طبق algorithm به مرحله ۲ برمی‌گردیم
   - ✅ اگر همگرا شد → گام فعلی تمام است — به سراغ گام بار بعدی (مرحله ۱) می‌رویم، تا رسیدن به تعداد گام‌های `analyze N`.

> 🚨 اگر تا سقف iteration همگرا نشود (مثلاً بار خیلی بزرگ برداشته شده یا مدل واقعاً ناپایدار شده)، `analyze` مقدار غیرصفر برمی‌گرداند و OpenSees به گام بعد نمی‌رود — باید در اسکریپت این را چک کرد (یا با گام بار کوچک‌تر، یا با تغییر algorithm، دوباره امتحان کرد).

### نقش algorithm در این چرخه

algorithm تعیین می‌کند مرحله‌ی ۲ (تشکیل tangent) چقدر دقیق و چقدر پرهزینه انجام شود:

| Algorithm | رفتار | سرعت همگرایی | هزینه هر iteration |
|---|---|---|---|
| `Newton` | tangent را در *هر* iteration دوباره می‌سازد | سریع (درجه ۲) | بالا |
| `ModifiedNewton` | فقط یک‌بار در ابتدای گام tangent می‌سازد و همان را در iteration‌های بعدی استفاده می‌کند | کندتر، iteration بیشتری لازم دارد | پایین |
| `KrylovNewton` | Newton را با شتاب‌دهنده‌ی Krylov ترکیب می‌کند | معمولاً سریع‌تر از ModifiedNewton | متوسط |

> 💡 **ارتباط با راهنمای build:** حالا معنای هشدار `-ffast-math` در مرحله‌ی ۴ راهنمای build روشن‌تر می‌شود: چون همین چرخه به مقایسه‌ی دقیق اعداد (residual کوچک‌تر از تلورانس؟) وابسته است، هر انحراف کوچک در محاسبات اعشاری می‌تواند مسیر همگرایی را عوض کند. به همین دلیل build با فلگ‌های IEEE-safe (`-O3` بدون `-ffast-math`) انجام شد.

---

*OpenSees 3.8.0 — Ubuntu 22.04 — Intel MKL 2026.0*
*Build تأیید‌شده روی Intel Core i9-13900K — راهنما شامل مفاهیم build و مکانیزم تحلیل خطی/غیرخطی*