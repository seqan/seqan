#!/usr/bin/env python2
"""Some validation for proc_doc.Proc*"""

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>'

class ProcDocValidator(object):
    """Validate proc_doc.Proc* objects.

    Implements the visitor pattern.
    """

    def __init__(self, msg_printer):
        self.msg_printer = msg_printer

    def validate(self, proc_entry):
        return


class MissingSignatureValidator(ProcDocValidator):
    """Validates for missing or empty signature."""

    def validate(self, proc_entry):
        IGNORED = ['variable', 'member_variable', 'tag', 'grouped_tag', 'typedef',
                   'grouped_typedef', 'signature', 'concept', 'member_typedef',
                   'enum', 'grouped_enum', 'enum_value']
        if not hasattr(proc_entry, 'signatures') or proc_entry.kind in IGNORED:
            return  # Skip if type has no signatures.
        if not proc_entry.signatures:
            msg = 'Missing @signature for this entry!'
            self.msg_printer.printTokenError(proc_entry.raw.first_token, msg, 'warning')


class MissingSignatureKeywordsValidator(ProcDocValidator):
    """Validates for missing keywords in signature (e.g. "class" for @class)."""

    def validate(self, proc_entry):
        if proc_entry.kind not in ['class', 'specialization']:
            return  # only handle those
        for i, sig in enumerate(proc_entry.raw.signatures):
            # TODO(holtgrew): Really allow typedef and ::Type/mfns here?
            if 'class ' not in sig.text.text and 'struct ' not in sig.text.text \
                    and 'typedef ' not in sig.text.text and 'using' not in sig.text.text \
                    and not sig.text.text.strip().endswith('::Type') \
                    and not sig.text.text.strip().endswith('::Type;'):
                msg = 'Missing keyword "class", "struct", "typedef", "using" in signature.'
                self.msg_printer.printTokenError(proc_entry.raw.signatures[i].text.tokens[0], msg, 'warning')


class OnlyRemarksInBodyValidator(ProcDocValidator):
    """Validates for the body starting with '@section Remarks'."""

    def validate(self, proc_entry):
        if not hasattr(proc_entry, 'body') or not proc_entry.body.children:
            return  # only handle if has non-empty body
        if proc_entry.body.children[0].type in ['h1', 'h2', 'h3', 'h4', 'h5'] and \
                proc_entry.body.children[0].children and \
                proc_entry.body.children[0].children[0].text == 'Remarks':
                msg = 'Detailed descrition starts with Remarks'
                self.msg_printer.printTokenError(proc_entry.raw.first_token, msg, 'warning')


class MissingParameterDescriptionValidator(ProcDocValidator):
    """Warns if the description is missing for a @param or @return."""

    def validate(self, proc_entry):
        if not hasattr(proc_entry, 'params') and \
           not hasattr(proc_entry, 'tparams') and \
           not hasattr(proc_entry, 'returns'):
            return  # Skip if type has no parameters
        # Check for empty name.
        for key in ['params', 'tparams', 'returns']:
            if not hasattr(proc_entry, key):
                continue  # Skip if missing.
            for val in getattr(proc_entry, key):
                if hasattr(val, 'name') and not val.name:
                    msg = 'Missing name for @%s' % key[:-1]
                elif hasattr(val, 'type') and not val.type:
                    msg = 'Missing type for @%s' % key[:-1]
                else:
                    continue  # skip
                self.msg_printer.printTokenError(val.raw.first_token, msg, 'warning')
        # Check for empty description.
        for key in ['params', 'tparams', 'returns']:
            if not hasattr(proc_entry, key):
                continue  # Skip if missing.
            for val in getattr(proc_entry, key):
                if val.desc.empty:
                    msg = 'Missing description for @%s' % key[:-1]
                    self.msg_printer.printTokenError(val.raw.first_token, msg, 'warning')


class ReturnVoidValidator(ProcDocValidator):
    """Warns if there is a (superfluous) @return void entry."""

    def validate(self, proc_entry):
        if not hasattr(proc_entry, 'returns'):
            return  # Skip if type has no returns member.
        for r in proc_entry.returns:
            if r.type == 'void':
                msg = '@return superfluous for "void" type -- simply show "void" in signature.'
                self.msg_printer.printTokenError(r.raw.first_token, msg, 'warning')


class EmptyBriefValidator(ProcDocValidator):
    """Warns if there is no non-empty @brief section for an entry."""

    def validate(self, proc_entry):
        IGNORED = ['mainpage', 'page']
        if proc_entry.kind in IGNORED:
            return  # Skip.
        if not hasattr(proc_entry, 'brief'):
            return  # Skip if type has no returns member.
        if not proc_entry.brief or proc_entry.brief.empty:
            msg = 'Missing non-empty @brief clause.'
            self.msg_printer.printTokenError(proc_entry.raw.first_token, msg, 'warning')


class ClassCannotExtendConceptValidator(ProcDocValidator):
    """Warns a class extends a concept."""

    def validate(self, proc_entry):
        if proc_entry.kind not in ['class', 'specialization']:
            return  # Skip.
        for ext in proc_entry.extends:
            if not ext in proc_entry.doc.top_level_entries:
                continue
            if proc_entry.doc.top_level_entries[ext].kind not in ['class', 'specialization']:
                msg = 'Class %s tries to inherit from non-class %s' % (proc_entry.name, ext)
                self.msg_printer.printTokenError(proc_entry.raw.first_token, msg, 'error')


# Array with the validator classes to use.
VALIDATORS = [MissingSignatureValidator,
              MissingParameterDescriptionValidator,
              MissingSignatureKeywordsValidator,
              #OnlyRemarksInBodyValidator,
              ReturnVoidValidator,
              EmptyBriefValidator,
              ClassCannotExtendConceptValidator]
