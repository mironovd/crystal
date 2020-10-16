# sage_import.py
#import imp
import types
import inspect
import os
import sys
import importlib.abc, importlib.util

import sage.all

def extendpath(path):
    sys.path.insert(0,os.path.abspath(path))

def extendpath_rel(path):
    extendpath(os.path.dirname(os.path.realpath(__file__))+"/"+path)

class StringLoader(importlib.abc.SourceLoader):
    def __init__(self, data, path):
        self.data = data
        self.path=path

    def get_source(self, fullname):
        return self.data

    def get_source(self, fullname):
        return self.data
    
    def get_data(self, path):
        return self.data.encode("utf-8")
    
    def get_filename(self, fullname):
        return self.path

def sage_import(modname, fromlist=None, namespace=None, newmodname=None):
    """
    Import a .sage module from the filename <modname>.sage

    Returns the resulting Python module.  If ``fromlist`` is given, returns
    just those members of the module into the global namespace where the
    function was called, or the given namespace.
    """

    filename = modname + '.sage'

    for path in sys.path:
        modpath = os.path.join(path, filename)
        if os.path.isfile(modpath):
            break
    else:
        raise ImportError('no file {} on sys.path'.format(filename))

    if newmodname is not None:
        modname = newmodname

    with open(modpath) as fobj:
        code = sage.all.preparse(fobj.read())
#        mod = imp.new_module(modname)
        loader = StringLoader(code,modpath)
        spec = importlib.util.spec_from_loader(modname, loader, origin=modpath)
        mod = importlib.util.module_from_spec(spec)
        sys.modules[modname] = mod
#        spec.loader.exec_module(mod)

#        mod = types.ModuleType(modname)
#        mod.__file__ = modpath
        # Fill with all the default Sage globals
        # We could just do a dict.update but we want to exclude dunder
        # and private attributes I guess
        for k, v in sage.all.__dict__.items():
            if not k.startswith('_'):
                mod.__dict__[k] = v

#        exec(code , mod.__dict__)
        spec.loader.exec_module(mod)

    if namespace is None:
        namespace = inspect.currentframe().f_back.f_globals


    if fromlist is not None:
        # First check that each name in fromlist exists before adding
        # any of them to the given namespace.
        for name in fromlist:
            if name not in mod.__dict__:
                raise ImportError('cannot import name {!r} from {}'.format(
                     name, filename))

        for name in fromlist:
            namespace[name] = mod.__dict__[name]
    else:
        namespace[modname] = mod

