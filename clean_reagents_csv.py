"""Clean non-standard Unicode characters from reagents.csv"""

import csv

file_path = r'c:\Git-softwares\Condition-agent\data\reagent_db\reagents.csv'
output_path = r'c:\Git-softwares\Condition-agent\data\reagent_db\reagents_cleaned.csv'

# Define character replacements
replacements = {
    # Corrupted dashes/hyphens (CJK characters that should be dashes)
    '鈥': '-',
    '慴': 'b',
    '慶': 'c',
    '慹': 'e',
    '慸': '-',
    '揹': '-',
    '揗': '-',
    '测': '-',
    '€': '-',
    '慍': 'C',
    
    # Middle dot variants
    '路': '·',
    '、': ',',
    
    # Chinese characters (likely corrupted chemical names)
    '乙': 'eth',
    '甲': 'meth',
    '苯': 'benz',
    '酮': 'one',
    '酸': 'acid',
    '氯': 'chloro',
    '溴': 'bromo',
    '碘': 'iodo',
    
    # Special spaces and marks
    '聽': ' ',
    '\x8e': '',
    
    # Smart quotes and primes
    '\u2019': "'",  # Right single quotation mark
    '\u2032': "'",  # Prime symbol,
    
    # Special symbols
    '«': '"',
    '®': '(R)',
    '±': '+/-',
    '³': '3',
    '÷': '/',
    
    # Latin extended characters (keep as-is or replace)
    'Í': 'I',
    'Ñ': 'N',
    'Ò': 'O',
    'à': 'a',
    'á': 'a',
    'æ': 'ae',
    'ñ': 'n',
    
    # Greek letters (lowercase - keep if valid scientific notation)
    'η': 'eta',
    'κ': 'kappa',
    'λ': 'lambda',
    'μ': 'mu',
    'Μ': 'M',
    
    # Replacement character (data corruption) - likely meant to be apostrophe/prime
    '\ufffd': "'",
}

# Read and process the file
with open(file_path, 'r', encoding='utf-8') as infile:
    content = infile.read()

original_content = content

# Apply replacements
for old_char, new_char in replacements.items():
    count = content.count(old_char)
    if count > 0:
        print(f"Replacing {repr(old_char)} with {repr(new_char)}: {count} occurrences")
    content = content.replace(old_char, new_char)

# Write the cleaned file
with open(output_path, 'w', encoding='utf-8', newline='') as outfile:
    outfile.write(content)

print(f'\nCleaned file saved to: {output_path}')

# Verify the changes
non_ascii_count = sum(1 for char in content if ord(char) > 127)
# Keep valid scientific Greek letters
valid_greek = set('αβγδεζηθικλμνξοπρστυφχψωΑΒΓΔΕΖΗΘΙΚΛΜΝΞΟΠΡΣΤΥΦΧΨΩ·')
problematic_chars = set(char for char in content if ord(char) > 127 and char not in valid_greek)
print(f'\nRemaining non-ASCII characters: {non_ascii_count}')
if problematic_chars:
    print(f'Problematic characters still present ({len(problematic_chars)}): {sorted(problematic_chars)}')
else:
    print('All problematic characters cleaned (only valid Greek letters and middle dots remain)')

# Show sample of changes
print("\n--- Sample Changes ---")
lines_changed = 0
for i, (orig_line, new_line) in enumerate(zip(original_content.split('\n'), content.split('\n')), 1):
    if orig_line != new_line:
        lines_changed += 1
        if lines_changed <= 5:
            print(f"\nLine {i}:")
            print(f"  BEFORE: {orig_line[:120]}")
            print(f"  AFTER:  {new_line[:120]}")

print(f"\nTotal lines changed: {lines_changed}")
