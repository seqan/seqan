#!/usr/bin/env python2
"""Code for translating a DDDoc tree and its node into raw_doc objects.
"""

import copy
import lexer
import raw_doc
import re
import sys


formatter = raw_doc.DoxFormatter(120)

class TokenTranslator(object):
    def translate(self, token_list):
        result = self.translateLinks(token_list)
        result = self.translateTT(result)
        return result

    def translateLinks(self, token_list):
        result = []
        for token in token_list:
            if '@' in token.val:
                vals = re.split('(@[^@]*@)', token.val)
                for val in vals:
                    if val.startswith('@') and val.endswith('@'):
                        t1 = copy.deepcopy(token)
                        t1.type = 'COMMAND_LINK'
                        t1.val = '@link'
                        t2 = copy.deepcopy(token)
                        if '|' in val:
                            type_name, title = val[1:-1].split('|')
                            type_name = translateTypename(type_name)
                            t2.val = '%s %s' % (type_name, title)
                        else:
                            t2.val = translateTypename(val[1:-1])
                        t2.val = ' %s ' % t2.val
                        t3 = copy.deepcopy(token)
                        t3.type = 'COMMAND_ENDLINK'
                        t3.val = '@endlink'
                        result += [t1, t2, t3]
                    else:
                        t = copy.deepcopy(token)
                        t.val = val
                        result.append(t)
            else:
                result.append(token)
        return result
        
    def translateTT(self, token_list):
        result = []
        for token in token_list:
            if '$' in token.val:
                vals = re.split('(\$[^\$]*\$)', token.val)
                for val in vals:
                    if val.startswith('$') and val.endswith('$'):
                        t1 = copy.deepcopy(token)
                        t1.type = 'HTML_TAG'
                        t1.val = '<tt>'
                        t2 = copy.deepcopy(token)
                        t2.val = val[1:-1]
                        t3 = copy.deepcopy(token)
                        t3.type = 'HTML_TAG'
                        t3.val = '</tt>'
                        result += [t1, t2, t3]
                    else:
                        t = copy.deepcopy(token)
                        t.val = val
                        result.append(t)
            else:
                result.append(token)
        return result

    
def translateTokens(tokens):
    return TokenTranslator().translate(tokens)


def translateTypename(name):
    if '.' in name:
        return name.split('.', 1)[1]
    else:
        return name


def translate(text):
    return text.replace('\colon', ':').replace('"', '')

# TODO(holtgrew): Translate spaces to underscores.


def migratePages(node):
    #print >>sys.stderr, 'Migrating pages...'
    pages = []
    #print node
    for name, child in node.children.iteritems():
        if name == 'Glossary':
            continue
        page = raw_doc.RawPage()
        page.title.tokens.append(lexer.Token('WORD', name, 0, 0, 0))
        s = 'Page'
        page.name.tokens.append(lexer.Token('WORD', s + page.title.text, 0, 0, 0))
        if child.children.get('summary'):
            for text in child.children['summary'].texts:
                if text.startswith('type=text:'):
                    t = lexer.Token('WORD', text[len('type=text:'):], 0, 0, 0)
                else:
                    t = lexer.Token('WORD', text, 0, 0, 0)
                raw_text = raw_doc.RawText(translateTokens([t]))
                raw_brief = raw_doc.RawBrief(raw_text)
                page.briefs.append(raw_brief)
        if child.children.get('description'):
            for text in child.children['description'].texts:
                if text.startswith('type=text:'):
                    t = lexer.Token('WORD', text[len('type=text:'):], 0, 0, 0)
                    raw_par = raw_doc.RawParagraph(raw_doc.RawText(translateTokens([t])))
                elif text.startswith('type=section:#'):
                    t = lexer.Token('WORD', text[len('type=section:#'):], 0, 0, 0)
                    raw_par = raw_doc.RawSection(raw_doc.RawText(translateTokens([t])))
                elif text.startswith('type=subsection:#.#'):
                    t = lexer.Token('WORD', text[len('type=subsection:#.#'):], 0, 0, 0)
                    raw_par = raw_doc.RawSection(raw_doc.RawText(translateTokens([t])), 1)
                else:
                    t = lexer.Token('WORD', text, 0, 0, 0)
                    raw_par = raw_doc.RawParagraph(raw_doc.RawText(translateTokens([t])))
                page.body.addParagraph(raw_par)
        pages.append(page)
    print 'RESULTING PAGES %s' % [p.title.text for p in pages]
    return pages


class GenericMigration(object):
    def __init__(self, node):
        self.node = node
        self.entry_class = None
        self.sep = '#'
        self.is_type = False

    def work(self):
        raw_entries = []
        #print self.node.key
        #if self.node.key == 'Memfunc':
        #    print 'GOT', self.node.key
        for name, child in self.node.children.iteritems():
            #if 'StringSet' in name or self.node.key == 'Memfunc':
            #    print 'GOT', self.node.key, name
            #print name
            #if name == 'Pattern':
            #    import pdb; pdb.set_trace()
            if self.node.key == 'Memfunc' or self.node.key == 'Typedef':
                name = name.replace('#', '::')
            name = translate(name)
            if self.entry_class is raw_doc.RawGroup and not child.children.get('tag') and not child.children.get('value'):
                entry = raw_doc.RawTag()
            else:
                entry = self.entry_class()
            if self.sep != '#' and '#' in name and not self.is_type:
                name = name.replace('#', self.sep)
            entry.name = raw_doc.RawText([lexer.Token('WORD', name, 0, 0, 0)])
            if self.entry_class is raw_doc.RawVariable:
                t = lexer.Token('WORD', 'VariableType', 0, 0, 0)
                entry.type = raw_doc.RawText([t])
            #if hasattr(entry, 'deprecation_msgs'):
            #    msg = raw_doc.RawDeprecated(raw_doc.RawText([lexer.Token('WORD', 'Fake deprecation message!', 0, 0, 0)]))
            #    entry.deprecation_msgs.append(msg)
            #print child
            #print name, child.children.keys()
            processed = set()
            if child.children.get('summary'):
                self.processSummary(entry, child.children['summary'])
                processed.add('summary')
            if child.children.get('general'):
                self.processExtends(entry, child.children['general'])
                processed.add('general')
            if child.children.get('implements'):
                self.processExtends(entry, child.children['implements'])
                processed.add('implements')
            if child.children.get('conceptimplements'):
#                self.processExtends(entry, child.children['conceptimplements'])
                processed.add('conceptimplements')
            if child.children.get('baseconcept'):
                self.processExtends(entry, child.children['baseconcept'])
                processed.add('baseconcept')
            if child.children.get('derived'):
#                self.processExtends(entry, child.children['derived'])
                processed.add('derived')
            if child.children.get('base'):
                self.processExtends(entry, child.children['base'])
                processed.add('base')
            if child.children.get('description'):
                self.processDescription(entry, child.children['description'])
                processed.add('description')
            if child.children.get('text'):
                self.processDescription(entry, child.children['text'])
                processed.add('text')
            if child.children.get('remarks'):
                self.processRemarks(entry, child.children['remarks'])
                processed.add('remarks')
            if child.children.get('remark'):
                self.processRemarks(entry, child.children['remark'], False)
                processed.add('remark')
            if child.children.get('note'):
                self.processRemarks(entry, child.children['note'])
                processed.add('note')
            if child.children.get('notes'):
                self.processRemarks(entry, child.children['notes'])
                processed.add('notes')
            first_example = True
            if child.children.get('example'):
                self.processExample(entry, child.children['example'])
                first_example = False
                processed.add('example')
            if child.children.get('file'):
                self.processFiles(entry, child.children['file'])
                processed.add('file')
            processed.add('cat')  # TODO(holtgrew): Handle category.
            processed.add('category')  # TODO(holtgrew): Handle category.
            if child.children.get('signature'):
                #print child
                self.processSignatures(entry, child.children['signature'])
                processed.add('signature')
            if child.children.get('demo'):
                self.processDemo(entry, child.children['demo'])
                processed.add('demo')
            if child.children.get('see'):
                self.processSee(entry, child.children['see'])
                processed.add('see')
            if child.children.get('returns'):
                self.processReturns(entry, child.children['returns'])
                processed.add('returns')
            if child.children.get('return'):
                self.processReturns(entry, child.children['return'])
                processed.add('return')
            if child.children.get('param'):
                self.processParams(entry, child.children['param'])
                processed.add('param')
            if child.children.get('params'):
                self.processParams(entry, child.children['params'])
                processed.add('params')
            if child.children.get('wiki'):
                self.processWiki(entry, child.children['wiki'])
                processed.add('wiki')
            if child.children.get('cite'):
                self.processCite(entry, child.children['cite'])
                processed.add('cite')
            if child.children.get('status'):
                self.processStatus(entry, child.children['status'])
                processed.add('status')
            if self.entry_class is raw_doc.RawEnum:
                if child.children.get('value'):
                    values = self.processValues(entry, child, child.children['value'])
                    #for v in values:
                    #    print v.getFormatted(formatter)
                    raw_entries += values
                    processed.add('value')
            if self.entry_class is raw_doc.RawGroup:  # tag
                if child.children.get('tag'):
                    values = self.processTags(entry, child, child.children['tag'])
                    #for v in values:
                    #    print v.getFormatted(formatter)
                    raw_entries += values
                    processed.add('tag')
                if child.children.get('value'):
                    values = self.processTags(entry, child, child.children['value'])
                    #for v in values:
                    #    print v.getFormatted(formatter)
                    raw_entries += values
                    processed.add('value')
            if child.children.get('include'):
                if hasattr(entry, 'addHeaderfile'):
                    self.processInclude(entry, child.children['include'])
                processed.add('include')
            if child.children.get('header'):
                if hasattr(entry, 'addHeaderfile'):
                    self.processInclude(entry, child.children['header'])
                processed.add('header')
            processed.add('classSpec')  # TODO(holtgrew): What's this?
            processed.add('default')  # TODO(holtgrew): What's this?
            processed.add('nowarn')  # TODO(holtgrew): What's this?
            processed.add('DISABLED')
            processed.add('hidefromindex')
            # We don't have double linking anymore.
            processed.add('function')
            processed.add('interfacefunc')
            processed.add('conceptfunc')
            processed.add('shortcutfor')
            processed.add('metafunction')
            processed.add('type')
            processed.add('spec')
            processed.add('memvar')
            processed.add('memfunc')
            processed.add('typedef')
            processed.add('shortcut')
            processed.add('metaraw_entries')
            processed.add('conceptmetafunc')
            processed.add('childconcept')
            processed.add('conceptusedbyfunc')
            if (child.children.get('class') or child.children.get('concept')) and not '#' in name and not '::' in name and not self.is_type and \
               self.entry_class in [raw_doc.RawFunction, raw_doc.RawMetafunction, raw_doc.RawVariable]:
                # Register raw_entries once for each class.
                texts = []
                if child.children.get('class'):
                    texts += child.children.get('class').texts
                if child.children.get('concept'):
                    texts += child.children.get('concept').texts
                texts = set(texts)
                for c in texts:
                    if '.' in c:
                        cname = translate(c.split('.', 1)[1])
                    else:
                        cname = translate(c)
                    entry2 = copy.deepcopy(entry)
                    entry2.name = raw_doc.RawText([lexer.Token('WORD', cname + self.sep + name, 0, 0, 0)])
                    raw_entries.append(entry2)
                    #print entry2.getFormatted(formatter)
            else:
                # Register raw_entries only once.
                raw_entries.append(entry)
                #print 'APPENDING', entry.getFormatted(formatter)
            processed.add('class')
            processed.add('concept')  # TODO(holtgrew): Ignoring for now...
            # Check that we processed all attributes.
            unhandled = set(child.children.keys()) - processed
            if unhandled:
                print 'Missed %s in %s (processed: %s)' % (unhandled, child, processed)
                sys.exit(1)
        return raw_entries

    def processSummary(self, entry, node):
        for text in node.texts:
            if text.startswith('type=text:'):
                t = lexer.Token('WORD', text[len('type=text:'):], 0, 0, 0)
                raw_text = raw_doc.RawText(translateTokens([t]))
            else:
                t = lexer.Token('WORD', text, 0, 0, 0)
                raw_text = raw_doc.RawText(translateTokens([t]))
            entry.addBrief(raw_doc.RawBrief(raw_text))

    def processExtends(self, entry, node):
        for text in node.texts:
            if '\u0001' in text:
                continue  # do not add inherited
            if text.startswith('Class.'):
                t = lexer.Token('WORD', text[len('Class.'):], 0, 0, 0)
                raw_text = raw_doc.RawText([t])
                entry.addExtends(raw_doc.RawExtends(raw_text))
            elif text.startswith('Spec.'):
                t = lexer.Token('WORD', text[len('Spec.'):], 0, 0, 0)
                raw_text = raw_doc.RawText([t])
                entry.addExtends(raw_doc.RawExtends(raw_text))
            elif text.startswith('Concept.'):
                t = lexer.Token('WORD', text[len('Concept.'):], 0, 0, 0)
                raw_text = raw_doc.RawText([t])
                if hasattr(entry, 'addImplements'):
                    entry.addImplements(raw_doc.RawImplements(raw_text))
                else:
                    entry.addExtends(raw_doc.RawExtends(raw_text))
            else:
                assert False, str(node)

    def processExample(self, entry, node):
        # Add example header.
        raw_text = raw_doc.RawText([lexer.Token('WORD', 'Examples', 0, 0, 0)])
        raw_section = raw_doc.RawSection(raw_text)
        entry.addParagraph(raw_section)
        for text in node.texts:
            if text.startswith('type=text:'):
                t = lexer.Token('WORD', text[len('type=text:'):], 0, 0, 0)
                raw_text = raw_doc.RawText(translateTokens([t]))
                raw_par = raw_doc.RawParagraph(raw_text)
            elif text.startswith('type=code:'):
                t = lexer.Token('WORD', '{.cpp}\n' + text[len('type=code:'):] + '\n', 0, 0, 0)
                raw_text = raw_doc.RawText(translateTokens([t]))
                raw_par = raw_doc.RawCode(raw_text)
            else:
                t = lexer.Token('WORD', text, 0, 0, 0)
                raw_text = raw_doc.RawText(translateTokens([t]))
                raw_par = raw_doc.RawParagraph(raw_text)
            entry.addParagraph(raw_par)

    def processFiles(self, entry, node):
        for text in node.texts:
            if text.startswith('../'):
                text = text[len('../'):]
            t = lexer.Token('WORD', text, 0, 0, 0)
            raw_text = raw_doc.RawText(translateTokens([t]))
            entry.addParagraph(raw_doc.RawInclude(raw_text))

    def processSee(self, entry, node):
        for text in node.texts:
            t = lexer.Token('WORD', translateTypename(text), 0, 0, 0)
            raw_text = raw_doc.RawText(translateTokens([t]))
            entry.addSee(raw_doc.RawSee(raw_text))

    def processValues(self, entry, node, value_nodes):
        result = []
        enum_name = translate(node.key)
        #print value_nodes
        for key, value in value_nodes.children.iteritems():
            raw_var = raw_doc.RawVariable()
            #print value
            var_name = value.key
            if '::' in enum_name:
                var_name = enum_name.split('::')[0] + '::' + var_name
            t = lexer.Token('WORD', enum_name, 0, 0, 0)
            raw_var.type = raw_doc.RawText([t])
            t = lexer.Token('WORD', enum_name + '#' + var_name, 0, 0, 0)
            raw_var.name = raw_doc.RawText([t])

            processed = set()
            if value.children.get('summary'):
                self.processSummary(raw_var, value.children['summary'])
                processed.add('summary')

            unhandled = set(value.children.keys()) - processed
            if unhandled:
                print 'Missed %s in %s' % (unhandled, child)
                sys.exit(1)

            #print raw_var.getFormatted(formatter)
            result.append(raw_var)
        return result

    def processTags(self, entry, node, tag_nodes):
        result = []
        group_name = translate(node.key)
        for key, value in tag_nodes.children.iteritems():
            raw_tag = raw_doc.RawTag()
            tag_name = value.key
            t = lexer.Token('WORD', group_name + '#' + tag_name, 0, 0, 0)
            raw_tag.name = raw_doc.RawText([t])

            processed = set()
            if value.children.get('summary'):
                self.processSummary(raw_tag, value.children['summary'])
                processed.add('summary')
            if value.children.get('remarks'):
                self.processRemarks(raw_tag, value.children['remarks'])
                processed.add('remarks')
            if value.children.get('remark'):
                self.processRemarks(raw_tag, value.children['remark'])
                processed.add('remark')
            if value.children.get('text'):
                self.processRemarks(raw_tag, value.children['text'])
                processed.add('text')
            if value.children.get('type'):
                ts = []
                if value.children['type'].texts:
                    ts.append(lexer.Token('WORD', ' Types: ', 0, 0, 0))
                    ts.append(lexer.Token('WORD', value.children['type'].texts[0], 0, 0, 0))
                for txt in value.children['type'].texts[1:]:
                    ts.append(lexer.Token('WORD', ', ' + txt, 0, 0, 0))
                raw_text = raw_doc.RawText(ts)
                entry.addParagraph(raw_doc.RawParagraph(raw_text))
                processed.add('type')

            processed.add('signature')  # Ignore.
            processed.add('include')  # TODO(holtgrew): Required here?
            processed.add('see')  # TODO(holtgrew): Required here?
            # We do not have double-linking any more.
            processed.add('function')

            unhandled = set(value.children.keys()) - processed
            if unhandled:
                print 'Missed %s in %s processed: %s' % (unhandled, node, processed)
                sys.exit(1)

            #print raw_tag.getFormatted(formatter)
            result.append(raw_tag)
        return result

    def processDemo(self, entry, node):
        for text in node.texts:
            t1 = lexer.Token('WORD', 'Demo: ', 0, 0, 0)
            t2 = lexer.Token('WORD', text, 0, 0, 0)
            raw_text = raw_doc.RawText(translateTokens([t1, t2]))
            entry.addParagraph(raw_doc.RawParagraph(raw_text))

    def processInclude(self, entry, node):
        for text in node.texts:
            t = lexer.Token('WORD', text, 0, 0, 0)
            raw_text = raw_doc.RawText(translateTokens([t]))
            raw_headerfile = raw_doc.RawHeaderfile(raw_text)
            entry.addHeaderfile(raw_headerfile)

    def processRemarks(self, entry, node, is_first=True):
        if is_first:
            raw_text = raw_doc.RawText([lexer.Token('WORD', 'Remarks', 0, 0, 0)])
            raw_section = raw_doc.RawSection(raw_text)
            entry.addParagraph(raw_section)
        for text in node.texts:
            if text.startswith('type=text:'):
                t = lexer.Token('WORD', text[len('type=text:'):], 0, 0, 0)
                raw_text = raw_doc.RawText(translateTokens([t]))
            else:
                t = lexer.Token('WORD', text, 0, 0, 0)
                raw_text = raw_doc.RawText(translateTokens([t]))
            entry.addParagraph(raw_doc.RawParagraph(raw_text))

    def processDescription(self, entry, node):
        for text in node.texts:
            if text.startswith('type=text:'):
                t = lexer.Token('WORD', text[len('type=text:'):], 0, 0, 0)
                raw_text = raw_doc.RawText(translateTokens([t]))
            else:
                t = lexer.Token('WORD', text, 0, 0, 0)
                raw_text = raw_doc.RawText(translateTokens([t]))
            entry.addParagraph(raw_doc.RawParagraph(raw_text))

    def processWiki(self, entry, node):
        for text in node.texts:
            t = lexer.Token('WORD', 'http://trac.seqan.de/wiki/' + text, 0, 0, 0)
            raw_text = raw_doc.RawText([t])
            entry.addParagraph(raw_doc.RawParagraph(raw_text))

    def processCite(self, entry, node):
        for text in node.texts:
            t = lexer.Token('WORD', text, 0, 0, 0)
            raw_text = raw_doc.RawText(translateTokens([t]))
            entry.addParagraph(raw_doc.RawParagraph(raw_text))

    def processStatus(self, entry, node):
        # TODO(holtgrew): Add support for @deprecated.
        for text in node.texts:
            t = lexer.Token('WORD', 'Status: ' + text, 0, 0, 0)
            raw_text = raw_doc.RawText(translateTokens([t]))
            entry.addParagraph(raw_doc.RawParagraph(raw_text))

    def processSignatures(self, entry, node):
        for text in node.texts:
            t = lexer.Token('WORD', text, 0, 0, 0)
            raw_text = raw_doc.RawText(translateTokens([t]))
            raw_sig = raw_doc.RawSignature(raw_text)
            entry.addSignature(raw_sig)

    def processParams(self, entry, node):
        for name, child in node.children.iteritems():
            ts = []
            if child.children.get('summary'):
                for summary in child.children['summary'].texts:
                    ts.append(lexer.Token('WORD', translate(summary), 0, 0, 0))
            if child.children.get('text'):
                for remark in child.children['text'].texts:
                    ts.append(lexer.Token('WORD', remark, 0, 0, 0))
            if child.children.get('tableheader'):
                for line in child.children['tableheader'].texts:
                    ts.append(lexer.Token('WORD', line, 0, 0, 0))
            if child.children.get('table'):
                for line in child.children['table'].texts:
                    ts.append(lexer.Token('WORD', line, 0, 0, 0))
            if child.children.get('remarks'):
                for remark in child.children['remarks'].texts:
                    ts.append(lexer.Token('WORD', remark, 0, 0, 0))
            if child.children.get('remark'):
                for remark in child.children['remark'].texts:
                    ts.append(lexer.Token('WORD', remark, 0, 0, 0))
            if child.children.get('node'):
                for node in child.children['node'].texts:
                    ts.append(lexer.Token('WORD', node, 0, 0, 0))
            if child.children.get('type'):
                if child.children['type'].texts:
                    ts.append(lexer.Token('WORD', ' Types: ', 0, 0, 0))
                    # TODO(holtgrew): Add @link?
                    ts.append(lexer.Token('WORD', translateTypename(child.children['type'].texts[0]), 0, 0, 0))
                for txt in child.children['type'].texts[1:]:
                    # TODO(holtgrew): Add @link?
                    ts.append(lexer.Token('WORD', ', ' + translateTypename(txt), 0, 0, 0))
            if child.children.get('concept'):
                if child.children['concept'].texts:
                    ts.append(lexer.Token('WORD', ' Concepts: ', 0, 0, 0))
                    ts.append(lexer.Token('WORD', child.children['concept'].texts[0], 0, 0, 0))
                for txt in child.children['concept'].texts[1:]:
                    ts.append(lexer.Token('WORD', ', ' + txt, 0, 0, 0))
            if child.children.get('class'):
                if child.children['class'].texts:
                    ts.append(lexer.Token('WORD', ' Classes: ', 0, 0, 0))
                    ts.append(lexer.Token('WORD', child.children['class'].texts[0], 0, 0, 0))
                for txt in child.children['class'].texts[1:]:
                    ts.append(lexer.Token('WORD', ', ' + txt, 0, 0, 0))
            if child.children.get('default'):
                if child.children['default'].texts:
                    ts.append(lexer.Token('WORD', ' Default: ', 0, 0, 0))
                    ts.append(lexer.Token('WORD', child.children['default'].texts[0], 0, 0, 0))
                for txt in child.children['default'].texts[1:]:
                    ts.append(lexer.Token('WORD', ', ' + txt, 0, 0, 0))
            if child.children.get('value'):
                if child.children['value'].texts:
                    ts.append(lexer.Token('WORD', ' Values: ', 0, 0, 0))
                    ts.append(lexer.Token('WORD', child.children['value'].texts[0], 0, 0, 0))
                for txt in child.children['value'].texts[1:]:
                    ts.append(lexer.Token('WORD', ', ' + txt, 0, 0, 0))
            if node.children.get('metafunction'):
                if node.children['metafunction'].texts:
                    ts.append(lexer.Token('WORD', ' Metafunctions: ', 0, 0, 0))
                    ts.append(lexer.Token('WORD', node.children['metafunction'].texts[0], 0, 0, 0))
                for txt in node.children['metafunction'].texts[1:]:
                    ts.append(lexer.Token('WORD', ', ' + txt, 0, 0, 0))
            if child.children.get('see'):
                if child.children['see'].texts:
                    ts.append(lexer.Token('WORD', ' Sees: ', 0, 0, 0))
                    ts.append(lexer.Token('WORD', child.children['see'].texts[0], 0, 0, 0))
                for txt in child.children['see'].texts[1:]:
                    ts.append(lexer.Token('WORD', ', ' + txt, 0, 0, 0))
            ts = translateTokens(ts)
            raw_text = raw_doc.RawText(ts)
            name_text = raw_doc.RawText([lexer.Token('WORD', name, 0, 0, 0)])
            if hasattr(entry, 'addParam'):
                raw_param = raw_doc.RawParam(name_text, raw_text)
                entry.addParam(raw_param)
            else:
                raw_param = raw_doc.RawTParam(name_text, raw_text)
                entry.addTParam(raw_param)

            # Check that we processed all attributes.
            unhandled = set(child.children.keys()) - set(['summary', 'remarks', 'type', 'default', 'text', 'concept', 'metafunction', 'remark', 'note', 'see', 'class', 'value', 'nowarn', 'tableheader', 'table'])
            if unhandled:
                print 'Missed %s in %s' % (unhandled, node)
                sys.exit(1)

    def processReturns(self, entry, node):
        ts = []
        if node.children.get('summary'):
            for text in node.children['summary'].texts:
                t = lexer.Token('WORD', text, 0, 0, 0)
                ts.append(t)
        if node.children.get('text'):
            for text in node.children['text'].texts:
                t = lexer.Token('WORD', text, 0, 0, 0)
                ts.append(t)
        if node.children.get('remarks'):
            for text in node.children['remarks'].texts:
                t = lexer.Token('WORD', text, 0, 0, 0)
                ts.append(t)
        if node.children.get('note'):
            for text in node.children['note'].texts:
                t = lexer.Token('WORD', text, 0, 0, 0)
                ts.append(t)
        if node.children.get('type'):
            if node.children['type'].texts:
                ts.append(lexer.Token('WORD', ' Types: ', 0, 0, 0))
                ts.append(lexer.Token('WORD', translateTypename(node.children['type'].texts[0]), 0, 0, 0))
            for txt in node.children['type'].texts[1:]:
                ts.append(lexer.Token('WORD', ', ' + translateTypename(txt), 0, 0, 0))
        if node.children.get('metafunction'):
            if node.children['metafunction'].texts:
                ts.append(lexer.Token('WORD', ' Metafunctions: ', 0, 0, 0))
                ts.append(lexer.Token('WORD', node.children['metafunction'].texts[0], 0, 0, 0))
            for txt in node.children['metafunction'].texts[1:]:
                ts.append(lexer.Token('WORD', ', ' + txt, 0, 0, 0))
        if node.children.get('default'):
            if node.children['default'].texts:
                ts.append(lexer.Token('WORD', ' Default: ', 0, 0, 0))
                ts.append(lexer.Token('WORD', node.children['default'].texts[0], 0, 0, 0))
            for txt in node.children['default'].texts[1:]:
                ts.append(lexer.Token('WORD', ', ' + txt, 0, 0, 0))
        ts = translateTokens(ts)
        name = raw_doc.RawText([lexer.Token('WORD', 'TReturn', 0, 0, 0)])
        raw_text = raw_doc.RawText(ts)
        raw_return = raw_doc.RawReturn(name, raw_text)
        entry.addReturn(raw_return)
        # Check that we processed all attributes.
        unhandled = set(node.children.keys()) - set(['summary', 'remarks', 'type', 'metafunction', 'text', 'param', 'note', 'default'])
        if unhandled:
            print 'Missed %s in %s' % (unhandled, node)
            sys.exit(1)
            
    def _processTextNode(self, node):
        texts = []
        for text in node.texts:
            if text.startswith('type=text:'):
                t = lexer.Token('WORD', text[len('type=text:'):], 0, 0, 0)
                raw_text = raw_doc.RawParagraph(raw_doc.RawText([t]))
            elif text.startswith('type=section:#'):
                t = lexer.Token('WORD', text[len('type=section:#'):], 0, 0, 0)
                raw_text = raw_doc.RawSection(raw_doc.RawText([t]))
            elif text.startswith('type=subsection:#.#'):
                t = lexer.Token('WORD', text[len('type=subsection:#.#'):], 0, 0, 0)
                raw_text = raw_doc.RawSection(raw_doc.RawText([t]), 1)
            else:
                t = lexer.Token('WORD', text, 0, 0, 0)
                raw_par = raw_doc.RawParagraph(raw_doc.RawText([t]))
            pars.append(raw_par)
        return pars


class FunctionMigration(GenericMigration):
    def __init__(self, node):
        GenericMigration.__init__(self, node)
        self.entry_class = raw_doc.RawFunction


class ClassMigration(GenericMigration):
    def __init__(self, node):
        GenericMigration.__init__(self, node)
        self.entry_class = raw_doc.RawClass
        self.is_type = True

    
class EnumMigration(GenericMigration):
    def __init__(self, node):
        GenericMigration.__init__(self, node)
        self.entry_class = raw_doc.RawEnum
        self.is_type = True

    
class ConceptMigration(GenericMigration):
    def __init__(self, node):
        GenericMigration.__init__(self, node)
        self.entry_class = raw_doc.RawConcept
        self.is_type = True

    
class MacroMigration(GenericMigration):
    def __init__(self, node):
        GenericMigration.__init__(self, node)
        self.entry_class = raw_doc.RawMacro

    
class MetafunctionMigration(GenericMigration):
    def __init__(self, node):
        GenericMigration.__init__(self, node)
        self.entry_class = raw_doc.RawMetafunction

    
class TagMigration(GenericMigration):
    def __init__(self, node):
        GenericMigration.__init__(self, node)
        self.entry_class = raw_doc.RawGroup
        self.is_type = True


class MemvarMigration(GenericMigration):
    def __init__(self, node):
        GenericMigration.__init__(self, node)
        self.entry_class = raw_doc.RawVariable
        self.sep = '::'


class TypedefMigration(GenericMigration):
    def __init__(self, node):
        GenericMigration.__init__(self, node)
        self.entry_class = raw_doc.RawTypedef
        self.sep = '::'
        self.is_type = True


class AdaptionMigration(GenericMigration):
    def __init__(self, node):
        GenericMigration.__init__(self, node)
        self.entry_class = raw_doc.RawAdaption
        self.sep = '::'
        self.is_type = True


def addGroups(raw_entries):
    """Add implicitely given groups to raw_entries."""
    # Collect known names.
    names = set()
    for e in raw_entries:
        name = e.name.text
        if '#' not in name and '::' not in name:
            names.add(name)
    # Collect unknown names.
    unknown = set()
    for e in raw_entries:
        name = e.name.text
        if '#' in name and not '::' in name:
            prefix = name.split('#', 1)[0]
            if prefix not in names:
                unknown.add(prefix)
    # Add groups.
    for name in unknown:
        group = raw_doc.RawGroup()
        t = lexer.Token('WORD', name, 0, 0, 0)
        group.name = raw_doc.RawText(translateTokens([t]))
        raw_entries.append(group)


def migrate(tree):
    print 'Class.keys', tree.root.children['Class'].children.keys()
    res = []
    res += migratePages(tree.root.children['Page'])
    res += migratePages(tree.root.children['Indexpage'])
    res += FunctionMigration(tree.root.children['Function']).work()
    res += FunctionMigration(tree.root.children['Memfunc']).work()
    res += ClassMigration(tree.root.children['Class']).work()
    res += ClassMigration(tree.root.children['Spec']).work()
    res += ConceptMigration(tree.root.children['Concept']).work()
    res += EnumMigration(tree.root.children['Enum']).work()
    res += MacroMigration(tree.root.children['Macro']).work()
    res += MetafunctionMigration(tree.root.children['Metafunction']).work()
    res += TagMigration(tree.root.children['Tag']).work()
    res += MemvarMigration(tree.root.children['Memvar']).work()
    res += TypedefMigration(tree.root.children['Typedef']).work()
    res += AdaptionMigration(tree.root.children['Adaption']).work()
    res += TypedefMigration(tree.root.children['Shortcut']).work()
    addGroups(res)
    #res += migratePages(tree.root.children['Demo'])
    #for x in res:
    #    print x.getFormatted(formatter)
    return res

