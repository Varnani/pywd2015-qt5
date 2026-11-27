"""Golden outputs produced by the exact release binaries identified in provenance.json."""

import hashlib,json,shutil

from pathlib import Path

import pytest

from src.helpers.wd_utils.wd_io import LCIO,DCIO

from src.helpers.wd_utils.wd_containers import LCParameterContainer,DCParameterContainer

ROOT=Path(__file__).parent/'fixtures/wd2015-win64'





def test_native_fixture_hashes():

    for name,expected in json.loads((ROOT/'provenance.json').read_text())['fixture_sha256'].items():

        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==expected





def test_native_lc_velocity_reference(tmp_path):

    io=LCIO(LCParameterContainer(),str(tmp_path),'unused')

    shutil.copyfile(ROOT/'lcout.active',io._get_output_path())

    table=io.read_synthetic_velocity_curve()

    assert len(table)==9 and len(table[0])==21

    assert table[1][:3]==[0,.05,.1]

    assert table[6][:3]==pytest.approx([-11.1978,-28.6190,-40.0949])





def test_native_dc_reference(tmp_path):

    params=DCParameterContainer();params['jdphs']=2

    io=DCIO(params,str(tmp_path),'unused')

    shutil.copyfile(ROOT/'dcout.active',io._get_output_path())

    result=io.read_results()

    assert result[0][4]==pytest.approx(5.0000011602)

    assert result[1][4]==pytest.approx(.4999956211)

    observations=io.read_unweighted_observations()

    assert len(observations)==6 and len(observations[0])==21

    assert io.read_solution_stats()[0][0]>0





@pytest.mark.parametrize('engine_type',['lc','dc'])

def test_real_engine_with_unicode_spaced_executable(app,tmp_path,engine_type):

    import os

    engine_dir=os.environ.get('PYWD_ENGINE_DIR')

    if not engine_dir:

        pytest.skip('Set PYWD_ENGINE_DIR for optional native binary integration checks.')

    directory=tmp_path/'motor Türkçe boşluk'

    directory.mkdir()

    for source in Path(engine_dir).iterdir():

        if source.suffix in ('.dat','.tab'):

            shutil.copyfile(source,directory/source.name)

    original=Path(engine_dir)/(engine_type+'_static_win64bit.exe')

    executable=directory/('hesap Türkçe '+engine_type+'.exe')

    shutil.copyfile(original,executable)

    params=LCParameterContainer() if engine_type=='lc' else DCParameterContainer()

    if engine_type=='dc':params['jdphs']=2

    io=(LCIO if engine_type=='lc' else DCIO)(params,str(directory),executable.name)

    io._input=(ROOT/(engine_type+'in.active')).read_text()

    io.timeout=30

    io.save().run()

    assert not (directory/(engine_type+'in.active')).exists()

    if engine_type=='lc':

        assert len(io.read_synthetic_velocity_curve()[0])==21

    else:

        assert io.read_results()[1][4]==pytest.approx(.4999956211,abs=1e-8)



@pytest.mark.parametrize('mode,ifoc,columns',[(1,1,7),(2,0,4),(1,0,5)])
def test_native_two_iteration_contract(tmp_path,mode,ifoc,columns):
    params=DCParameterContainer();params['jdphs']=mode;params['ifoc']=ifoc;params.keeps['niter']=2
    io=DCIO(params,str(tmp_path),'unused')
    fixture=ROOT/f'dc-mode{mode}-ifoc{ifoc}-niter2'/'dcout.active'
    shutil.copyfile(fixture,io._get_output_path())
    assert len(io.read_unweighted_observations())==columns
    assert len(io.read_unweighted_observations()[0])==21
    assert io.read_results()[1][4]==pytest.approx(.5,abs=1e-5)


def test_truncated_final_solution_never_borrows_previous_dimensions(tmp_path):
    params=DCParameterContainer();params['jdphs']=2;params.keeps['niter']=2
    io=DCIO(params,str(tmp_path),'unused')
    text=(ROOT/'dcout.active').read_text()
    second=text[:text.index('  2   pole')]
    Path(io._get_output_path()).write_text(text+second)
    with pytest.raises(ValueError):io.read_component_dimensions()


def test_truncated_final_solution_never_borrows_previous_stats(tmp_path):
    params=DCParameterContainer();params['jdphs']=2;params.keeps['niter']=2
    io=DCIO(params,str(tmp_path),'unused')
    text=(ROOT/'dcout.active').read_text()
    second=text[:text.index('   Mean residual for input values')]
    Path(io._get_output_path()).write_text(text+second)
    with pytest.raises(ValueError):io.read_solution_stats()


def test_curve_stats_match_final_native_iteration(tmp_path):
    params = DCParameterContainer()
    params.keeps['niter'] = 2
    io = DCIO(params, str(tmp_path), 'unused')
    fixture = ROOT / 'dc-mode1-ifoc1-niter2' / 'dcout.active'
    shutil.copyfile(fixture, io._get_output_path())
    stats = io.read_curve_stats()
    assert stats[0][1] == 21
    assert stats[0][2] == pytest.approx(0.000025119343)


@pytest.mark.parametrize('replacement', [None, '      1           21        invalid'])
def test_final_curve_stats_cannot_borrow_previous_or_publish_invalid(tmp_path, replacement):
    params = DCParameterContainer()
    params.keeps['niter'] = 2
    io = DCIO(params, str(tmp_path), 'unused')
    text = (ROOT / 'dc-mode1-ifoc1-niter2' / 'dcout.active').read_text()
    header = '    Curve   No. of obs.        Std. dev.'
    start = text.rfind(header)
    end = text.index('\n \n', start)
    if replacement is None:
        text = text[:start] + text[end:]
    else:
        row_start = text.index('\n', start) + 1
        text = text[:row_start] + replacement + text[end:]
    Path(io._get_output_path()).write_text(text)
    with pytest.raises(ValueError):
        io.read_curve_stats()


@pytest.mark.parametrize('mode,nparams,columns', [(2, 35, 39), (1, 37, 42)])
@pytest.mark.parametrize('execute', [False, True], ids=['golden', 'native'])
def test_native_wide_dc_records(tmp_path, mode, nparams, columns, execute):
    import os
    engine_dir = os.environ.get('PYWD_ENGINE_DIR')
    if execute and not engine_dir:
        pytest.skip('Set PYWD_ENGINE_DIR for optional native binary integration checks.')
    params = DCParameterContainer()
    params['jdphs'] = mode
    io = DCIO(params, engine_dir or str(tmp_path), 'dc_static_win64bit.exe')
    fixture = ROOT / f'dc-mode{mode}-ifoc1-wide{nparams}'
    if execute:
        io._input = (fixture / 'dcin.active').read_text()
        io.timeout = 180
        io.save().run()
    else:
        shutil.copyfile(fixture / 'dcout.active', io._get_output_path())
    results = io.read_results()
    assert len(results) == nparams
    assert [row[0] for row in results] == [56] * nparams
    assert [row[1] for row in results] == list(range(1, nparams + 1))
    assert [row[4] for row in results] == pytest.approx([1.0086008038] * nparams, abs=1e-8)
    table = io.read_unweighted_observations()
    assert len(table) == columns
    assert all(len(column) == 21 * nparams for column in table)
    assert table[0][:21] == pytest.approx([i * (0.1 if mode == 1 else 0.05) for i in range(21)])
    assert len(io.read_curve_stats()) == nparams
