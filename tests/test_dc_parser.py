from pathlib import Path
import pytest
from src.helpers.wd_utils.wd_io import DCIO
from src.helpers.wd_utils.wd_containers import DCParameterContainer


def fixture_io(tmp_path, nparams=35, jdphs=2, ifoc=1, width=38, truncate=False):
    params = DCParameterContainer()
    params["jdphs"] = jdphs
    params["ifoc"] = ifoc
    params.keeps["niter"] = 1
    io = DCIO(params, str(tmp_path), "unused")
    base = (2 if jdphs == 1 else 1) + 2 * ifoc + 1
    count = nparams + base
    result = "Input-Output in F Format\n\nNo. Curve Input Correction Output Sigma\n"
    for n in range(nparams):
        result += f"{n + 1:5d}{0:4d}{1.0:19.10f}{0.1:18.10f}{1.1:19.10f}{0.01:18.10f}\n"
    equation = "Unweighted Observational Equations\n\nheader\n"
    for row in range(2):
        values = [float(i + 1 + row * 100) for i in range(count)]
        if truncate and row == 1:
            values.pop()
        for start in range(0, len(values), width):
            chunk = values[start : start + width]
            equation += (
                f"{chunk[0]:13.5f}" + "".join(f"{v:20.10f}" for v in chunk[1:])
                if jdphs == 1
                else "".join(f"{v:20.10f}" for v in chunk)
            ) + "\n"
    Path(io._get_output_path()).write_text(equation + "\n" + result + "\n")
    return io, count


@pytest.mark.parametrize(
    "nparams,jdphs,ifoc,width",
    [
        (34, 2, 1, 38),
        (35, 2, 1, 38),
        (36, 1, 1, 41),
        (37, 1, 1, 41),
        (35, 2, 0, 38),
        (37, 1, 0, 41),
        (35, 2, 1, 20),
        (37, 1, 1, 23),
        (59, 2, 1, 20),
    ],
)
def test_logical_equations_follow_parameter_count(tmp_path, nparams, jdphs, ifoc, width):
    io, count = fixture_io(tmp_path, nparams, jdphs, ifoc, width)
    table = io.read_unweighted_observations()
    assert len(table) == count
    assert table[0] == [1, 101]
    assert table[-1] == [count, count + 100]


def test_truncated_equations_rejected(tmp_path):
    io, _ = fixture_io(tmp_path, truncate=True)
    with pytest.raises(ValueError):
        io.read_unweighted_observations()


def test_equations_match_same_iteration_result_block(tmp_path):
    io, count = fixture_io(tmp_path)
    old = Path(io._get_output_path()).read_text()
    newer, newcount = fixture_io(tmp_path, nparams=3)
    Path(io._get_output_path()).write_text(old + Path(newer._get_output_path()).read_text())
    io.parameters.keeps["niter"] = 2
    assert len(io.read_unweighted_observations()) == newcount


def test_dimensions_from_final_solution(tmp_path):
    params = DCParameterContainer()
    io = DCIO(params, str(tmp_path), "unused")

    def row(star, value):
        return f"{star:3d}{'pole':>7}{value:14.6f}{0.1:14.6f}{0.2:14.6f}{0.01:14.6f}\n"

    Path(io._get_output_path()).write_text(
        "Input-Output in F Format\n"
        + row(1, 1)
        + "\n"
        + row(2, 2)
        + "\nInput-Output in F Format\n"
        + row(1, 3)
        + "\n"
        + row(2, 4)
        + "\n"
    )
    dimensions = io.read_component_dimensions()
    assert dimensions[0][2] == [3]
    assert dimensions[1][2] == [4]
