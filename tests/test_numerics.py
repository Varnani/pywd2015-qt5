from datetime import date
import math
import pytest
from src.helpers import methods
from src.helpers.wd_utils.wd_containers import _ParameterContainer

@pytest.mark.parametrize("seed", [138472370, 100000000, 138472375, -100000000, 10, -10])
def test_integer_precision_preserves_significant_zeroes(seed):
    field = _ParameterContainer.Parameter("seed", float, seed).format(11, 0, "F")
    assert len(field) == 11
    assert float(field) == seed

@pytest.mark.parametrize("value,precision", [(1e-10, 5), (1.2e-20, 5), (10.125, 5), (-0.125, 5), (0, 0)])
def test_fractional_and_exponent_fields(value, precision):
    field = _ParameterContainer.Parameter("p", float, value).format(20, precision, "D")
    assert len(field) == 20
    assert float(field.replace("D", "e")) == pytest.approx(value)

@pytest.mark.parametrize("ymd", [(1900,3,1),(1901,1,1),(2100,1,1),(2000,1,1),(2024,2,29),(1582,10,15)])
def test_gregorian_jd_against_independent_date_delta(ymd):
    expected = 2440587.5 + (date(*ymd)-date(1970,1,1)).days
    assert methods.convert_ut_to_jd(*ymd,0,0,0) == expected
    assert methods.convert_jd_to_ut(expected)[:3] == ymd

def test_j2000_noon():
    assert methods.convert_ut_to_jd(2000,1,1,12,0,0) == 2451545.0

@pytest.mark.parametrize("start,end", [(.15,1),(-1.5,-.5),(0,0),(0,3.2),(-.8,.3)])
def test_alias_preserves_all_observations(start,end):
    x,y = [.1,.2,.9],[1,2,9]
    expected=sorted((v+k,w) for v,w in zip(x,y) for k in range(-5,6) if start <= v+k <= end)
    for actual in [methods.alias_phased_obs_with_phase(x,y,start,end), methods.alias_phased_obs_with_jd(x,y,start,end,0,1)]:
        assert sorted(zip(*actual)) == expected

@pytest.mark.parametrize("p", [0,-1,float("nan"),float("inf")])
def test_invalid_period_rejected(p):
    with pytest.raises(ValueError):
        methods.alias_phased_obs_with_jd([],[],0,1,0,p)

def test_empty_alias():
    assert methods.alias_phased_obs_with_phase([],[],0,1) == ([],[])


def test_historical_julian_leap_day_roundtrip():
    jd=methods.convert_ut_to_jd(1500,2,29,0,0,0)
    assert methods.convert_jd_to_ut(jd)[:3]==(1500,2,29)


def test_alias_exact_endpoints_against_bruteforce():
    import random
    rng=random.Random(20261125)
    for _ in range(100):
        phases=[rng.random() for _ in range(3)]
        period=rng.uniform(.1,3)
        t0=2451545.0
        base=[t0+x*period for x in phases]
        start=base[0]-period
        stop=base[2]+2*period
        if stop<start: start,stop=stop,start
        expected=sorted((x+k*period,index) for index,x in enumerate(base) for k in range(-5,6) if start<=x+k*period<=stop)
        actual=methods.alias_phased_obs_with_jd(phases,[0,1,2],start,stop,t0,period)
        assert sorted(zip(*actual))==expected


@pytest.mark.parametrize('value',[1e-6,-1e-6,.00001234])
def test_fortran_fixed_fields_never_emit_an_f_exponent(value):
    field=_ParameterContainer.Parameter('fixed',float,value).format(12,6,'F')
    assert 'F' not in field
    assert float(field)==pytest.approx(round(value,6),abs=1e-12)
