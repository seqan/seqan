#!/usr/env/bin python
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
                   'enum', 'grouped_enum']
        if not hasattr(proc_entry, 'signatures') or proc_entry.kind in IGNORED:
            return  # Skip if type has no signatures.
        if not proc_entry.signatures:
            msg = 'Missing @signature for this entry!'
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
    """Warns if there is a (superflous) @return void entry."""

    def validate(self, proc_entry):
        if not hasattr(proc_entry, 'returns'):
            return  # Skip if type has no returns member.
        for r in proc_entry.returns:
            if r.type == 'void':
                msg = '@return superflous for "void" type -- simply show "void" in signature.'
                self.msg_printer.printTokenError(r.raw.first_token, msg, 'warning')

# Array with the validator classes to use.
VALIDATORS = [MissingSignatureValidator,
              MissingParameterDescriptionValidator,
              ReturnVoidValidator]
