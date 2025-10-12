# Windows Environment Variable Setup Guide

## Setting API Keys Permanently in Windows

### Quick Method (PowerShell - Recommended)

```powershell
# Set OpenAI API key (permanent for current user)
[System.Environment]::SetEnvironmentVariable('OPENAI_API_KEY', 'sk-your-key-here', 'User')

# Set Aliyun/DeepSeek API key (optional)
[System.Environment]::SetEnvironmentVariable('ALIYUN_API_KEY', 'sk-your-key-here', 'User')

# IMPORTANT: Close and reopen PowerShell after setting!

# Verify in new PowerShell window:
echo $env:OPENAI_API_KEY
echo $env:ALIYUN_API_KEY
```

### GUI Method (Visual)

1. **Open System Properties**
   - Press `Win + R`
   - Type: `sysdm.cpl`
   - Press Enter

2. **Navigate to Environment Variables**
   - Click **"Advanced"** tab
   - Click **"Environment Variables"** button at bottom

3. **Add New Variable**
   - Under **"User variables for [YourName]"**
   - Click **"New..."** button
   - Variable name: `OPENAI_API_KEY`
   - Variable value: `sk-your-actual-key-here`
   - Click **OK**

4. **Save Changes**
   - Click **OK** on Environment Variables window
   - Click **OK** on System Properties window

5. **Restart Terminal**
   - Close all PowerShell/CMD windows
   - Open a new one to load the new variable

### Command Prompt Method (setx)

```cmd
REM Set the variable
setx OPENAI_API_KEY "sk-your-key-here"

REM You'll see: SUCCESS: Specified value was saved.

REM Close and reopen CMD/PowerShell!
```

---

## Testing Your Setup

### 1. Check Variable is Set

```powershell
# In a NEW PowerShell window
echo $env:OPENAI_API_KEY

# Should output: sk-...
# If it shows nothing, the variable isn't set or you're in the same session
```

### 2. Test with Python

```powershell
# Quick Python test
python -c "import os; print('API Key:', os.getenv('OPENAI_API_KEY', 'NOT FOUND'))"

# Should show: API Key: sk-...
```

### 3. Run LLMTools Examples

```powershell
cd C:\Git-softwares\Condition-agent
python llmtools/examples.py

# Should show:
# ✓ OPENAI_API_KEY found: sk-xxxxxxxx...
```

---

## Common Issues

### Issue 1: "Variable not found" in new terminal

**Cause**: Still in old PowerShell session

**Fix**: 
```powershell
# Close PowerShell completely
# Open NEW PowerShell window
echo $env:OPENAI_API_KEY
```

### Issue 2: Works in PowerShell but not Python

**Cause**: Python started before variable was set

**Fix**:
```powershell
# Restart Python/IDE
# Or reload in current Python session:
import os
os.environ['OPENAI_API_KEY'] = os.getenv('OPENAI_API_KEY')
```

### Issue 3: Variable disappears after reboot

**Cause**: Set with `$env:` instead of `SetEnvironmentVariable`

**Fix**: Use the permanent method above

---

## Temporary Setting (Session Only)

If you just want to test without permanent storage:

```powershell
# Current session only (lost when window closes)
$env:OPENAI_API_KEY = "sk-your-key-here"

# Verify immediately (same session)
echo $env:OPENAI_API_KEY

# Use right away
python llmtools/examples.py
```

---

## Security Best Practices

### 1. Never Commit API Keys

Add to `.gitignore`:
```gitignore
# API Keys
.env
.env.local
**/api_keys.txt
```

### 2. Use .env Files for Projects (Alternative)

Create `.env` file in project root:
```bash
# .env
OPENAI_API_KEY=sk-your-key-here
ALIYUN_API_KEY=sk-your-key-here
```

Install python-dotenv:
```powershell
pip install python-dotenv
```

The llmtools will automatically load it!

### 3. Rotate Keys Regularly

- Regenerate API keys periodically
- Update environment variables when keys change
- Monitor API usage for suspicious activity

---

## Quick Command Reference

```powershell
# Set permanent (User level)
[System.Environment]::SetEnvironmentVariable('VAR_NAME', 'value', 'User')

# Set permanent (System level - requires admin)
[System.Environment]::SetEnvironmentVariable('VAR_NAME', 'value', 'Machine')

# Get current value
[System.Environment]::GetEnvironmentVariable('VAR_NAME', 'User')

# Delete variable
[System.Environment]::SetEnvironmentVariable('VAR_NAME', $null, 'User')

# List all user variables
[System.Environment]::GetEnvironmentVariables('User')

# Reload variable in current session
$env:VAR_NAME = [System.Environment]::GetEnvironmentVariable('VAR_NAME', 'User')
```

---

## Next Steps

1. **Set your API key** using Method 1 or 2 above
2. **Close and reopen PowerShell**
3. **Verify** with `echo $env:OPENAI_API_KEY`
4. **Test** with `python llmtools/examples.py`
5. **Start using** llmtools in your workflows!

---

**Note**: If you're using Git Bash, WSL, or Linux subsystems, the commands are different. Let me know if you need those instructions.
