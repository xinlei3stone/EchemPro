#!/usr/bin/env bash
# 重建 EchemPro.app —— 源码正本位于本目录 (EchemPro-src/)
#
# 用法:
#   ./build.sh                 # 复用 RuHeterojunctionHER/.venv
#   PYTHON=/path/to/python ./build.sh
#
# 说明: 旧 app 会先改名为 ../EchemPro.app.<时间戳>.bak, 不直接删除。

set -euo pipefail
cd "$(dirname "$0")"

PY="${PYTHON:-../RuHeterojunctionHER/.venv/bin/python}"
if [ ! -x "$PY" ]; then
    echo "找不到可执行的 Python: $PY" >&2
    echo "请设置 PYTHON=/path/to/python (需已安装 requirements.txt 中的依赖)" >&2
    exit 1
fi

echo "==> PyInstaller ($("$PY" --version 2>&1))"
"$PY" -m PyInstaller --windowed --name EchemPro \
    --workpath build --noconfirm EchemPro.py

APP="dist/EchemPro.app"
if [ ! -d "$APP" ]; then
    echo "构建失败: 未生成 $APP" >&2
    exit 1
fi

if [ -d ../EchemPro.app ]; then
    BAK="../EchemPro.app.$(date +%Y%m%d-%H%M%S).bak"
    echo "==> 备份旧应用 -> $BAK"
    mv ../EchemPro.app "$BAK"
fi

mv "$APP" ../EchemPro.app

# BUNDLE 之外, PyInstaller 还会留下同内容的 onedir 目录 (dist/EchemPro/), 属于构建中间产物
if [ -d dist/EchemPro ]; then
    echo "==> 清理构建中间产物 dist/EchemPro"
    rm -rf dist/EchemPro
fi

echo "==> 完成: $(cd .. && pwd)/EchemPro.app"
