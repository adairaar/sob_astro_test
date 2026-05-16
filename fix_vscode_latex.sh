#!/bin/sh
# fix_vscode_latex.sh — wires VS Code LaTeX Workshop to the full TeX Live 2026

SETTINGS="$HOME/Library/Application Support/Code/User/settings.json"

echo "=== Searching for xelatex installations ==="
# Search for binaries in bin/ subdirectories (includes symlinks via -L)
FOUND=$(find -L /usr/local/texlive /Library/TeX /usr/texbin /usr/local/bin \
    -name "xelatex" 2>/dev/null \
    | grep -E "bin/[^/]+/xelatex|bin/xelatex" \
    | sort -r)

if [ -z "$FOUND" ]; then
    echo "ERROR: No xelatex found."
    exit 1
fi

echo "Found:"
echo "$FOUND" | sed 's/^/  /'

# Prefer the newest full (non-basic) installation
BEST=$(echo "$FOUND" | grep -v "basic" | grep "2026" | head -1)
[ -z "$BEST" ] && BEST=$(echo "$FOUND" | grep -v "basic" | grep "2024" | head -1)
[ -z "$BEST" ] && BEST=$(echo "$FOUND" | grep -v "basic" | head -1)
[ -z "$BEST" ] && BEST=$(echo "$FOUND" | head -1)

echo ""
echo "Selected: $BEST"

# Verify it's executable
if [ ! -x "$BEST" ]; then
    echo "ERROR: $BEST is not executable"
    exit 1
fi

# Check threeparttable availability
KPSEWHICH="$(dirname "$BEST")/kpsewhich"
if [ -x "$KPSEWHICH" ]; then
    TPT=$("$KPSEWHICH" threeparttable.sty 2>/dev/null)
    if [ -n "$TPT" ]; then
        echo "threeparttable.sty found: $TPT  ✓"
    else
        echo "threeparttable.sty: not in this TeX tree (project stub will be used)"
    fi
fi

# ── Patch VS Code settings via Python ─────────────────────────────────────
echo ""
echo "=== Patching VS Code settings ==="
mkdir -p "$(dirname "$SETTINGS")"
[ -f "$SETTINGS" ] || echo "{}" > "$SETTINGS"

python3 - "$SETTINGS" "$BEST" <<'PYEOF'
import sys, json, pathlib

settings_path = pathlib.Path(sys.argv[1])
xelatex = sys.argv[2]
bibtex   = str(pathlib.Path(xelatex).parent / "bibtex")

try:
    settings = json.loads(settings_path.read_text())
except Exception:
    settings = {}

settings["latex-workshop.latex.tools"] = [
    {
        "name": "xelatex",
        "command": xelatex,
        "args": ["-synctex=1", "-interaction=nonstopmode",
                 "-file-line-error", "%DOC%"]
    },
    {
        "name": "bibtex",
        "command": bibtex,
        "args": ["%DOCFILE%"]
    }
]
settings["latex-workshop.latex.recipes"] = [
    {"name": "XeLaTeX",          "tools": ["xelatex"]},
    {"name": "XeLaTeX + BibTeX", "tools": ["xelatex","bibtex","xelatex","xelatex"]}
]
settings["latex-workshop.latex.recipe.default"] = "XeLaTeX"

settings_path.write_text(json.dumps(settings, indent=2))
print("  Wrote: " + str(settings_path))
PYEOF

echo ""
echo "=== Done ==="
echo "  1. Quit VS Code completely  (Cmd+Q)"
echo "  2. Reopen VS Code"
echo "  3. Open sob_paper.tex and press Cmd+Alt+B to build"
