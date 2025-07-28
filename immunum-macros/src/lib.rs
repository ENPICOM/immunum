use proc_macro::TokenStream;
use quote::quote;
use syn::{parse_macro_input, Data, DeriveInput};

/// Custom derive macro that generates FromStr implementations with alias support
#[proc_macro_derive(ParseFromString, attributes(value))]
pub fn derive_parse_from_string(input: TokenStream) -> TokenStream {
    let input = parse_macro_input!(input as DeriveInput);
    let name = &input.ident;

    let expanded = match &input.data {
        Data::Enum(data_enum) => {
            let mut match_arms = Vec::new();

            for variant in &data_enum.variants {
                let variant_name = &variant.ident;
                let variant_str = variant_name.to_string().to_lowercase();

                // Start with the main variant name
                let mut patterns = vec![variant_str.clone()];

                // Extract aliases from clap's #[value(alias = "...")] attributes
                for attr in &variant.attrs {
                    if attr.path().is_ident("value") {
                        // Convert the attribute meta to string and parse manually
                        let meta_str = quote! { #attr }.to_string();

                        // Simple regex-like parsing for alias = "value" patterns
                        if meta_str.contains("alias") {
                            // Extract quoted strings after alias =
                            let parts: Vec<&str> = meta_str.split("alias").collect();
                            for part in parts.iter().skip(1) {
                                // Look for = "value" pattern
                                if let Some(eq_pos) = part.find('=') {
                                    let after_eq = &part[eq_pos + 1..].trim();
                                    if let Some(start) = after_eq.find('"') {
                                        if let Some(end) = after_eq[start + 1..].find('"') {
                                            let alias = &after_eq[start + 1..start + 1 + end];
                                            patterns.push(alias.to_lowercase());
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                // Create match arm for all patterns
                let pattern_strs: Vec<_> = patterns.iter().map(|p| quote! { #p }).collect();
                match_arms.push(quote! {
                    #(#pattern_strs)|* => Ok(#name::#variant_name),
                });
            }

            // Collect all valid values for error messages
            let mut valid_values = Vec::new();
            for variant in &data_enum.variants {
                let variant_name = &variant.ident;
                let variant_str = variant_name.to_string();
                let mut aliases = Vec::new();

                // Extract aliases from clap's #[value(alias = "...")] attributes
                for attr in &variant.attrs {
                    if attr.path().is_ident("value") {
                        let meta_str = quote! { #attr }.to_string();
                        if meta_str.contains("alias") {
                            let parts: Vec<&str> = meta_str.split("alias").collect();
                            for part in parts.iter().skip(1) {
                                if let Some(eq_pos) = part.find('=') {
                                    let after_eq = &part[eq_pos + 1..].trim();
                                    if let Some(start) = after_eq.find('"') {
                                        if let Some(end) = after_eq[start + 1..].find('"') {
                                            let alias = &after_eq[start + 1..start + 1 + end];
                                            aliases.push(alias.to_string());
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                // Format as "VARIANT (alias1, alias2)" or just "VARIANT" if no aliases
                let formatted_value = if aliases.is_empty() {
                    variant_str
                } else {
                    format!("{} ({})", variant_str, aliases.join(", "))
                };
                valid_values.push(formatted_value);
            }

            let valid_values_str = valid_values.join(", ");
            let type_name = name.to_string();

            // Generate the FromStr implementation
            quote! {
                impl std::str::FromStr for #name {
                    type Err = String;

                    fn from_str(s: &str) -> Result<Self, Self::Err> {
                        match s.to_lowercase().as_str() {
                            #(#match_arms)*
                            _ => Err(format!("{} not supported: '{}', use any of: {}",
                                #type_name, s, #valid_values_str)),
                        }
                    }
                }

                impl #name {
                    /// Parse from string with case-insensitive matching and alias support
                    pub fn parse_from_string(s: &str) -> Result<Self, String> {
                        s.parse()
                    }

                    /// Parse a vector of strings into a vector of this type
                    pub fn parse_vec_from_strings(strings: Vec<String>) -> Result<Vec<Self>, String> {
                        strings.into_iter()
                            .map(|s| Self::parse_from_string(&s))
                            .collect()
                    }
                }
            }
        }
        _ => {
            return syn::Error::new_spanned(
                &input,
                "ParseFromString can only be derived for enums",
            )
            .to_compile_error()
            .into();
        }
    };

    TokenStream::from(expanded)
}
