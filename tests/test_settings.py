from configparser import ConfigParser
from pathlib import Path
from src import constants
from src.interfaces.lcdcpicker_interface import Widget


def test_literal_paths_and_legacy_fallback(app,tmp_path,monkeypatch):
    legacy=tmp_path/'legacy-config';user=tmp_path/'user'/'pywd_config.dat'
    executable=tmp_path/'100% Türkçe binary.exe';executable.write_text('fixture')
    parser=ConfigParser(interpolation=None);parser['wd_paths']={'lc_path':str(executable),'dc_path':str(executable)}
    with legacy.open('w',encoding='utf-8') as out:parser.write(out)
    monkeypatch.setattr(constants,'PYWD_CONFIG_PATH',str(user))
    monkeypatch.setattr(constants,'CONFIG_DIR',str(user.parent))
    monkeypatch.setattr(constants,'LEGACY_CONFIG_PATH',str(legacy))
    widget=Widget()
    monkeypatch.setattr(widget,'exec_',lambda:(_ for _ in ()).throw(AssertionError('valid legacy config ignored')))
    assert widget.check_config() is True
    widget.save()
    reloaded=ConfigParser(interpolation=None);reloaded.read(user,encoding='utf-8')
    assert reloaded.get('wd_paths','lc_path')==str(executable)
