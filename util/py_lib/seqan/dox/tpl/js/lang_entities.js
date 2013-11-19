function darken(hex, amount) {
	var color = new less.tree.Color(hex.substr(1), 1.0);
	amount = new less.tree.Value(amount);
	return less.tree.functions.darken(color, amount).toCSS();
}

window.colors = {{ config.colors|tojson }};

window.langEntities = {{ config.lang_entities|tojson }};