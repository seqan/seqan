#!/usr/bin/env python2

import sys
import unittest

import seqan.dox.sig_parser as sig_parser


class TestParseEnum(unittest.TestCase):
    def testValid(self):
        txt = 'enum EnumName'
        parser = sig_parser.SigParser(txt)
        entry = parser.parse()
        self.assertEqual(entry.kind, 'enum')
        self.assertEqual(entry.name, 'EnumName')
        txt = 'enum EnumName;'
        self.assertEqual(entry.toString(), txt)


class TestParseStruct(unittest.TestCase):
    def testValid(self):
        txt = 'struct StructName'
        parser = sig_parser.SigParser(txt)
        entry = parser.parse()
        self.assertEqual(entry.kind, 'struct')
        self.assertEqual(entry.name, 'StructName')
        txt = 'struct StructName;'
        self.assertEqual(entry.toString(), txt)


class TestParseClass(unittest.TestCase):
    def testValid(self):
        txt = 'class ClassName'
        parser = sig_parser.SigParser(txt)
        entry = parser.parse()
        self.assertEqual(entry.kind, 'class')
        self.assertEqual(entry.name, 'ClassName')
        txt = 'class ClassName;'
        self.assertEqual(entry.toString(), txt)


class TestParseConcept(unittest.TestCase):
    def testValid(self):
        txt = 'concept ConceptConcept'
        parser = sig_parser.SigParser(txt)
        entry = parser.parse()
        self.assertEqual(entry.kind, 'concept')
        self.assertEqual(entry.name, 'ConceptConcept')
        txt = 'concept ConceptConcept;'
        self.assertEqual(entry.toString(), txt)


class TestParseVariable(unittest.TestCase):
    def testValid(self):
        txt = 'Type variable;'
        parser = sig_parser.SigParser(txt)
        entry = parser.parse()
        self.assertEqual(entry.kind, 'variable')
        self.assertEqual(entry.var_type, 'Type')
        self.assertEqual(entry.name, 'variable')
        txt = 'Type variable;'
        self.assertEqual(entry.toString(), txt)


class TestParseFunction(unittest.TestCase):
    def testEmptyParams(self):
        txt = 'void foo()'
        parser = sig_parser.SigParser(txt)
        entry = parser.parse()
        self.assertEqual(entry.kind, 'function')
        self.assertEqual(entry.return_type, 'void')
        self.assertEqual(entry.is_tpl, False)
        self.assertEqual(len(entry.params), 0)

    def testWithParams(self):
        txt = 'void foo(int x, double y)'
        parser = sig_parser.SigParser(txt)
        entry = parser.parse()
        self.assertEqual(entry.kind, 'function')
        self.assertEqual(entry.return_type, 'void')
        self.assertEqual(entry.is_tpl, False)
        self.assertEqual(len(entry.params), 2)
        self.assertEqual(entry.params[0].type, 'int')
        self.assertEqual(entry.params[0].name, 'x')
        self.assertEqual(entry.params[1].type, 'double')
        self.assertEqual(entry.params[1].name, 'y')
        txt = 'void foo(int x, double y);'
        self.assertEqual(entry.toString(), txt)

    def testMemberFunction(self):
        txt = 'void Foo::bar(int x, double y)'
        parser = sig_parser.SigParser(txt)
        entry = parser.parse()
        self.assertEqual(entry.kind, 'function')
        self.assertEqual(entry.name, 'Foo::bar')
        self.assertEqual(entry.return_type, 'void')
        self.assertEqual(entry.is_tpl, False)
        self.assertEqual(len(entry.params), 2)
        self.assertEqual(entry.params[0].type, 'int')
        self.assertEqual(entry.params[0].name, 'x')
        self.assertEqual(entry.params[1].type, 'double')
        self.assertEqual(entry.params[1].name, 'y')
        txt = 'void Foo::bar(int x, double y);'
        self.assertEqual(entry.toString(), txt)

    def testConstructorNoParams(self):
        txt = 'String::String()'
        parser = sig_parser.SigParser(txt)
        entry = parser.parse()
        self.assertEqual(entry.kind, 'function')
        self.assertEqual(entry.name, 'String::String')
        self.assertEqual(entry.return_type, None)
        self.assertEqual(entry.is_tpl, False)
        self.assertEqual(len(entry.params), 0)
        txt = 'String::String();'
        self.assertEqual(entry.toString(), txt)

    def testConstructorParams(self):
        txt = 'String::String(int x, double y)'
        parser = sig_parser.SigParser(txt)
        entry = parser.parse()
        self.assertEqual(entry.kind, 'function')
        self.assertEqual(entry.name, 'String::String')
        self.assertEqual(entry.return_type, None)
        self.assertEqual(entry.is_tpl, False)
        self.assertEqual(len(entry.params), 2)
        self.assertEqual(entry.params[0].type, 'int')
        self.assertEqual(entry.params[0].name, 'x')
        self.assertEqual(entry.params[1].type, 'double')
        self.assertEqual(entry.params[1].name, 'y')
        txt = 'String::String(int x, double y);'
        self.assertEqual(entry.toString(), txt)

    def testDestructorNoParams(self):
        txt = 'String::~String()'
        parser = sig_parser.SigParser(txt)
        entry = parser.parse()
        self.assertEqual(entry.kind, 'function')
        self.assertEqual(entry.name, 'String::~String')
        self.assertEqual(entry.return_type, None)
        self.assertEqual(entry.is_tpl, False)
        self.assertEqual(len(entry.params), 0)
        txt = 'String::~String();'
        self.assertEqual(entry.toString(), txt)

    def testInterfaceFunction(self):
        txt = 'void Foo#bar(int x, double y)'
        parser = sig_parser.SigParser(txt)
        entry = parser.parse()
        self.assertEqual(entry.kind, 'function')
        self.assertEqual(entry.name, 'Foo#bar')
        self.assertEqual(entry.return_type, 'void')
        self.assertEqual(entry.is_tpl, False)
        self.assertEqual(len(entry.params), 2)
        self.assertEqual(entry.params[0].type, 'int')
        self.assertEqual(entry.params[0].name, 'x')
        self.assertEqual(entry.params[1].type, 'double')
        self.assertEqual(entry.params[1].name, 'y')
        txt = 'void Foo#bar(int x, double y);'
        self.assertEqual(entry.toString(), txt)

    def testConstructor(self):
        pass


class TestTemplateFunction(unittest.TestCase):
    def testEmptyParams(self):
        txt = ('template <typename T1, int I>\n'
               'void foo()')
        parser = sig_parser.SigParser(txt)
        entry = parser.parse()
        self.assertEqual(entry.kind, 'function')
        self.assertEqual(entry.return_type, 'void')
        self.assertEqual(entry.is_tpl, True)
        self.assertEqual(len(entry.params), 0)
        self.assertEqual(len(entry.tparams), 2)
        self.assertEqual(entry.tparams[0].type, 'typename')
        self.assertEqual(entry.tparams[0].name, 'T1')
        self.assertEqual(entry.tparams[1].type, 'int')
        self.assertEqual(entry.tparams[1].name, 'I')
        txt = ('template <typename T1, int I>\n'
               'void foo();')
        self.assertEqual(entry.toString(), txt)

    def testWithParams(self):
        txt = ('template <typename T1, int I>\n'
               'void foo(int x, double y)')
        parser = sig_parser.SigParser(txt)
        entry = parser.parse()
        self.assertEqual(entry.kind, 'function')
        self.assertEqual(entry.return_type, 'void')
        self.assertEqual(entry.is_tpl, True)
        self.assertEqual(len(entry.params), 2)
        self.assertEqual(entry.params[0].type, 'int')
        self.assertEqual(entry.params[0].name, 'x')
        self.assertEqual(entry.params[1].type, 'double')
        self.assertEqual(entry.params[1].name, 'y')
        self.assertEqual(len(entry.tparams), 2)
        self.assertEqual(entry.tparams[0].type, 'typename')
        self.assertEqual(entry.tparams[0].name, 'T1')
        self.assertEqual(entry.tparams[1].type, 'int')
        self.assertEqual(entry.tparams[1].name, 'I')
        txt = ('template <typename T1, int I>\n'
               'void foo(int x, double y);')
        self.assertEqual(entry.toString(), txt)


class TestTemplateClass(unittest.TestCase):
    def testEmptyParams(self):
        txt = ('template <>\n'
               'class C')
        parser = sig_parser.SigParser(txt)
        entry = parser.parse()
        self.assertEqual(entry.kind, 'class')
        self.assertEqual(entry.name, 'C')
        self.assertEqual(entry.is_tpl, True)
        self.assertEqual(len(entry.tparams), 0)
        txt = ('template <>\n'
               'class C;')
        self.assertEqual(entry.toString(), txt)


    def testWithParams(self):
        txt = ('template <typename T1, int I>\n'
               'class C')
        parser = sig_parser.SigParser(txt)
        entry = parser.parse()
        self.assertEqual(entry.kind, 'class')
        self.assertEqual(entry.is_tpl, True)
        self.assertEqual(len(entry.tparams), 2)
        self.assertEqual(entry.name, 'C')
        self.assertEqual(entry.tparams[0].type, 'typename')
        self.assertEqual(entry.tparams[0].name, 'T1')
        self.assertEqual(entry.tparams[1].type, 'int')
        self.assertEqual(entry.tparams[1].name, 'I')
        txt = ('template <typename T1, int I>\n'
               'class C;')
        self.assertEqual(entry.toString(), txt)


class TestMetafunction(unittest.TestCase):
    def testValueNoParam(self):
        txt = 'TInt Metafunction<>::VALUE'
        parser = sig_parser.SigParser(txt)
        entry = parser.parse()
        self.assertEqual(entry.kind, 'metafunction')
        self.assertEqual(entry.name, 'Metafunction')
        self.assertEqual(entry.return_type, 'TInt')
        self.assertEqual(entry.return_name, 'VALUE')
        self.assertEqual(len(entry.tparams), 0)
        txt = 'TInt Metafunction<>::VALUE;'
        self.assertEqual(entry.toString(), txt)

    def testValueInterfaceNoParam(self):
        txt = 'TInt Klass#Metafunction<>::VALUE'
        parser = sig_parser.SigParser(txt)
        entry = parser.parse()
        self.assertEqual(entry.kind, 'metafunction')
        self.assertEqual(entry.name, 'Klass#Metafunction')
        self.assertEqual(entry.return_type, 'TInt')
        self.assertEqual(entry.return_name, 'VALUE')
        self.assertEqual(len(entry.tparams), 0)
        txt = 'TInt Klass#Metafunction<>::VALUE;'
        self.assertEqual(entry.toString(), txt)

    def testValueWithParam(self):
        txt = 'TInt Metafunction<T1, T2>::VALUE'
        parser = sig_parser.SigParser(txt)
        entry = parser.parse()
        self.assertEqual(entry.kind, 'metafunction')
        self.assertEqual(entry.name, 'Metafunction')
        self.assertEqual(entry.return_type, 'TInt')
        self.assertEqual(entry.return_name, 'VALUE')
        self.assertEqual(len(entry.tparams), 2)
        self.assertEqual(entry.tparams[0].name, 'T1')
        self.assertEqual(entry.tparams[1].name, 'T2')
        txt = 'TInt Metafunction<T1, T2>::VALUE;'
        self.assertEqual(entry.toString(), txt)

    def testTypeNoParam(self):
        txt = 'Metafunction<>::Type'
        parser = sig_parser.SigParser(txt)
        entry = parser.parse()
        self.assertEqual(entry.kind, 'metafunction')
        self.assertEqual(entry.name, 'Metafunction')
        self.assertEqual(entry.return_name, 'Type')
        self.assertEqual(len(entry.tparams), 0)
        txt = 'Metafunction<>::Type;'
        self.assertEqual(entry.toString(), txt)

    def testTypeInterfaceNoParam(self):
        txt = 'Klass#Metafunction<>::Type'
        parser = sig_parser.SigParser(txt)
        entry = parser.parse()
        self.assertEqual(entry.kind, 'metafunction')
        self.assertEqual(entry.name, 'Klass#Metafunction')
        self.assertEqual(entry.return_name, 'Type')
        self.assertEqual(len(entry.tparams), 0)
        txt = 'Klass#Metafunction<>::Type;'
        self.assertEqual(entry.toString(), txt)

    def testTypeWithParam(self):
        txt = 'Metafunction<T1, T2>::Type'
        parser = sig_parser.SigParser(txt)
        entry = parser.parse()
        self.assertEqual(entry.kind, 'metafunction')
        self.assertEqual(entry.name, 'Metafunction')
        self.assertEqual(entry.return_name, 'Type')
        self.assertEqual(len(entry.tparams), 2)
        self.assertEqual(entry.tparams[0].name, 'T1')
        self.assertEqual(entry.tparams[1].name, 'T2')
        txt = 'Metafunction<T1, T2>::Type;'
        self.assertEqual(entry.toString(), txt)


if __name__ == '__main__':
    unittest.main()
