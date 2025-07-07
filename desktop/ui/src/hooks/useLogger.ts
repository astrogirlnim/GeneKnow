import { useCallback } from 'react'

interface LogLevel {
  DEBUG: 'debug'
  INFO: 'info'
  WARN: 'warn'
  ERROR: 'error'
}

// Define a more specific type for log data
type LogData = Record<string, unknown> | unknown[] | string | number | boolean | null | undefined

// Removed unused LOG_LEVELS constant to fix TypeScript error

export const useLogger = () => {
  const log = useCallback((level: keyof LogLevel, message: string, data?: LogData) => {
    const timestamp = new Date().toISOString()
    const logMessage = `[${timestamp}] [${level}] ${message}`
    
    switch (level) {
      case 'DEBUG':
        console.debug(logMessage, data)
        break
      case 'INFO':
        console.info(logMessage, data)
        break
      case 'WARN':
        console.warn(logMessage, data)
        break
      case 'ERROR':
        console.error(logMessage, data)
        break
      default:
        console.log(logMessage, data)
    }
  }, [])

  const debug = useCallback((message: string, data?: LogData) => {
    log('DEBUG', message, data)
  }, [log])

  const info = useCallback((message: string, data?: LogData) => {
    log('INFO', message, data)
  }, [log])

  const warn = useCallback((message: string, data?: LogData) => {
    log('WARN', message, data)
  }, [log])

  const error = useCallback((message: string, data?: LogData) => {
    log('ERROR', message, data)
  }, [log])

  return {
    log,
    debug,
    info,
    warn,
    error
  }
} 