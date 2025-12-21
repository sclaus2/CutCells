"""Test CLI functionality."""

import subprocess
import sys
from pathlib import Path

# Add tools directory to path for imports
tools_dir = Path(__file__).parent.parent.parent
sys.path.insert(0, str(tools_dir))


def test_cli_list():
    """Test CLI list command."""
    result = subprocess.run(
        [sys.executable, "-m", "tools.cutcells_tables.cli", "list"],
        capture_output=True,
        text=True,
        cwd=tools_dir.parent,
    )

    assert result.returncode == 0
    assert "generate-hex" in result.stdout


def test_cli_dump_hex_base():
    """Test CLI dump-hex-base command."""
    result = subprocess.run(
        [sys.executable, "-m", "tools.cutcells_tables.cli", "dump-hex-base", "--mask", "173"],
        capture_output=True,
        text=True,
        cwd=tools_dir.parent,
    )

    assert result.returncode == 0
    assert "Mask: 173" in result.stdout
    assert "Base ID:" in result.stdout
